# Copyright (c) 2022: Oscar Dowson and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    convert_files(; directory::String)
"""
function convert_files(; directory::String)
    for (root, _, files) in walkdir(directory)
        println(root)
        for file in files
            if endswith(file, ".gms")
                println("  ", file)
                convert_file(joinpath(root, file))
            end
        end
    end
    return
end

"""
    convert_file(filename)

Convert a .gms file into the .nl equivalent.
"""
function convert_file(filename)
    model = _parse_gms_file(filename)
    nl_file = replace(filename, ".gms" => ".nl")
    MOI.write_to_file(model, nl_file)
    return
end

_replace_expr(x, expr) = expr

function _replace_expr(x::Vector{MOI.VariableIndex}, expr::Symbol)
    if expr == :POWER || expr == :power
        return :^
    elseif expr == :arctan
        return :atan
    elseif expr == :objvar
        return x[end]
    end
    m = match(r"x([0-9]+)", "$expr")
    if m === nothing
        return expr
    end
    return x[parse(Int, m[1])]
end

function _replace_expr(x, expr::Expr)
    if Meta.isexpr(expr, :call) && expr.args[1] == :sqr
        expr.args[1] = :^
        push!(expr.args, 2)
    end
    for i in 1:length(expr.args)
        expr.args[i] = _replace_expr(x, expr.args[i])
    end
    return expr
end

function _read_next_block(io)
    line = strip(readuntil(io, ';'))
    while startswith(line, "*")
        # Some comments are not `;` terminated. Strip them out and return the
        # following block.
        i = findfirst('\n', line)
        if i === nothing
            return line
        end
        line = string(strip(line[i+1:end]))

    end
    return line
end

function _parse_gms_file(filename)
    model = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    nlp = MOI.Nonlinear.Model()
    x = MOI.VariableIndex[]
    io = open(filename, "r")
    while !eof(io)
        line = _read_next_block(io)
        if startswith(line, '*')
            continue
        elseif startswith(line, "Variables")
            # variables are comma-separated.
            num_variables = 1 + count(isequal(','), line)
            for i in 1:num_variables
                xi = MOI.add_variable(model)
                push!(x, xi)
                name = i == num_variables ? "objvar" : "x$i"
                MOI.set(model, MOI.VariableName(), xi, name)
            end
        elseif startswith(line, "Positive Variables")
            for indices in findall(r"x([0-9]+)?", line)
                m = match(r"x([0-9]+)?", line[indices])
                xi = MOI.VariableIndex(parse(Int, m[1]))
                MOI.add_constraint(model, xi, MOI.GreaterThan(0.0))
            end
        elseif startswith(line, "Negative Variables")
            for indices in findall(r"x([0-9]+)?", line)
                m = match(r"x([0-9]+)?", line[indices])
                xi = MOI.VariableIndex(parse(Int, m[1]))
                MOI.add_constraint(model, xi, MOI.LessThan(0.0))
            end
        elseif startswith(line, "Equations")
            continue  # We don't need to pre-allocate variables.
        elseif startswith(line, r"e[0-9]+\.\.")
            # Isolate expression from name
            str_expr = String(split(line, ".. ")[2])
            # Remove new lines
            str_expr = replace(str_expr, "\n" => "")
            # Simplify expressions like + x - y + z into +(x, -y, z)
            str_expr = replace(str_expr, "- " => "+ -")
            # Correct (in)equalities.
            for pair in ["=L=" => "<=", "=G=" => ">=", "=E=" => "=="]
                str_expr = replace(str_expr, pair)
            end
            expr = Meta.parse(str_expr)
            _replace_expr(x, expr)
            func = expr.args[2]
            set = if expr.args[1] == :(<=)
                MOI.LessThan{Float64}(expr.args[3])
            elseif expr.args[1] == :(>=)
                MOI.GreaterThan{Float64}(expr.args[3])
            else
                @assert expr.args[1] == :(==)
                MOI.EqualTo{Float64}(expr.args[3])
            end
            try
                MOI.Nonlinear.add_constraint(nlp, func, set)
            catch err
                rethrow(err)
            end
        elseif startswith(line, "Solve") || startswith(line, "\$")
            if occursin("minimizing objvar", line)
                MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
                MOI.set(model, MOI.ObjectiveFunction{typeof(x[end])}(), x[end])
            elseif occursin("maximizing objvar", line)
                MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
                MOI.set(model, MOI.ObjectiveFunction{typeof(x[end])}(), x[end])
            else
                error("Unsupported objective: $line")
            end
        elseif isempty(line)
            continue
        elseif startswith(line, "Model")
            continue
        elseif startswith(line, "m.")
            continue
        elseif match(r"e([0-9]+).m", line) !== nothing
            continue  # Ignore non-default marginals
        else
            xi = if startswith(line, "objvar.")
                m = match(r"(objvar).([a-z]+) = (.+)", line)
                x[end]
            else
                m = match(r"x([0-9]+)\.([a-z]+) = (.+)", line)
                if m === nothing
                    error("unknown line: $line")
                end
                x[parse(Int, m[1])]
            end
            value = parse(Float64, m[3])
            if m[2] == "l"
                MOI.set(model, MOI.VariablePrimalStart(), xi, value)
            elseif m[2] == "lo"
                _set_or_modify(model, xi, MOI.GreaterThan(value))
            elseif m[2] == "up"
                _set_or_modify(model, xi, MOI.LessThan(value))
            elseif m[2] == "fx"
                _set_or_modify(model, xi, MOI.LessThan(value))
                _set_or_modify(model, xi, MOI.GreaterThan(value))
            else
                error(m[2])
            end
        end
    end
    evaluator =
        MOI.Nonlinear.Evaluator(nlp, MOI.Nonlinear.SparseReverseMode(), x)
    block = MOI.NLPBlockData(evaluator)
    MOI.set(model, MOI.NLPBlock(), block)
    return model
end

function _set_or_modify(model, f, s)
    ci = MOI.ConstraintIndex{typeof(f),typeof(s)}(f.value)
    if MOI.is_valid(model, ci)
        MOI.set(model, MOI.ConstraintSet(), ci, s)
    else
        MOI.add_constraint(model, f, s)
    end
    return
end

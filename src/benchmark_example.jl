# Copyright (c) 2022: Oscar Dowson and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

const SETUPS = Dict(
    "ASL" => (
        optimizer = () -> AmplNLWriter.Optimizer(Ipopt_jll.amplexe),
        backend = MOI.Nonlinear.SparseReverseMode(),
    ),
    "JuMP" => (
        optimizer = Ipopt.Optimizer,
        backend = MOI.Nonlinear.SparseReverseMode(),
    ),
    "SymbolicAD" => (
        optimizer = Ipopt.Optimizer,
        backend = MathOptSymbolicAD.DefaultBackend(),
    ),
)

"""
    benchmark_example(file::String, setup::String)
"""
function benchmark_example(file::String, setup::String)
    println("Benchmarking: $file")
    config = SETUPS[setup]
    model = JuMP.read_from_file(file)
    JuMP.set_optimizer(model, config.optimizer)
    JuMP.set_optimizer_attribute(model, "max_wall_time", 60.0)
    JuMP.set_optimizer_attribute(model, "print_level", 0)
    start = time()
    JuMP.optimize!(model; _differentiation_backend = config.backend)
    time_solve = time() - start
    status = JuMP.termination_status(model)
    objective = status == MOI.LOCALLY_SOLVED ? JuMP.objective_value(model) : NaN
    nlp_block = JuMP.MOI.get(model, JuMP.MOI.NLPBlock())
    nvar, ncon = _model_size(file)
    data = Dict(
        "file" => file,
        "backend" => setup,
        "variables" => nvar,
        "constraints" => ncon,
        "status" => status,
        "cost" => objective,
        "time_solve" => time_solve,
        "time_callbacks" =>
            nlp_block.evaluator.eval_objective_timer +
            nlp_block.evaluator.eval_objective_gradient_timer +
            nlp_block.evaluator.eval_constraint_timer +
            nlp_block.evaluator.eval_constraint_jacobian_timer +
            nlp_block.evaluator.eval_hessian_lagrangian_timer,
    )
    open("results.jsonl", "a") do io
        println(io, JSON.json(data))
    end
    return
end

function _model_size(file)
    open(file, "r") do io
        readline(io)
        line = split.(readline(io); keepempty = false)
        items = parse.(Int, String.(line))
        return items[1], items[2]
    end
end
# Copyright (c) 2022: Oscar Dowson and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

import AmplNLWriter
import CSV
import DataFrames
import Ipopt
import Ipopt_jll
import JSON
import JuMP
import MathOptInterface
import MathOptSymbolicAD
import Plots
import Statistics

const MOI = MathOptInterface

const SKIPS = Dict(
    "ASL" => String[],
    "JuMP" => String[],
    "SymbolicAD" => String[
        "gams/cute/arglina.nl",
        "gams/cute/chebyqad.nl",
        "gams/cute/watson.nl",
        "gams/cute/vanderm1.nl",
        "gams/cute/vanderm2.nl",
        "gams/cute/vanderm3.nl",
        "gams/jeffc/3a.nl",
        "gams/nnls/nnls.nl",
    ],
)

include("benchmark_example.jl")
include("convert_files.jl")
include("plots.jl")

function list_of_models(; directory::String, setup::String)
    list = String[]
    for (root, _, files) in walkdir(directory)
        for file in files
            if endswith(file, ".nl") && all(s -> !occursin(s, file), SKIPS[setup])
                push!(list, joinpath(root, file))
            end
        end
    end
    return list
end

const GAMS_DIR = joinpath(@__DIR__, "..", "gams")

if length(ARGS) > 0
    if ARGS[1] == "--convert"
        convert_files(directory = GAMS_DIR)
    elseif ARGS[1] == "--benchmark"
        setup = ARGS[2]
        for file in list_of_models(directory = GAMS_DIR, setup = setup)
            benchmark_example(file, setup)
        end
    end
end

# plot_comparisons(
#     "opf-solve.md",
#     [
#         "JuMP NL" => "ASL",
#         "CasADi" => "CasADi",
#         "OPOMO" => "Gravity",
#         "JuMP" => "JuMP",
#         "JuMP SymAD" => "SymbolicAD",
#     ],
# )
#
# plot_comparisons(
#     "opf-callbacks.md",
#     [
#         "JuMP" => "JuMP",
#         "JuMP SymAD" => "SymbolicAD",
#     ],
# )

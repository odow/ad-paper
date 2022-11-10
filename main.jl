import CSV
import DataFrames
import JSON
import Plots
import Statistics

function plot_solved(plot, data; name, label)
    labels = ["ASL", "CasADi", "Gravity", "JuMP", "SymbolicAD"]
    i = findfirst(isequal(label), labels)
    styles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    colors = Plots.cgrad(:darktest, 5, categorical = true)
    x = sort(data[!, name])
    y = (1:length(x)) ./ length(x)
    Plots.plot!(
        plot,
        x,
        y;
        label = label,
        linewidth = 2,
        linestyle = styles[i],
        color = colors[i],
    )
    return plot
end

function plot_comparisons(filename, fields)
    filename = joinpath(@__DIR__, filename)
    data = CSV.read(filename, DataFrames.DataFrame; skipto = 3)
    DataFrames.rename!(strip, data)
    DataFrames.select!(data, DataFrames.names(data)[2:end-1])
    plt = Plots.plot(;
        # xaxis = :log,
        legend = :bottomright,
        xlabel = "Time (s)",
        ylabel = "Proportion of problems solved",
        size = (450, 300),
    )
    for (name, label) in fields
        plot_solved(plt, data; name = name, label = label)
    end
    Plots.savefig(replace(filename, ".md" => ".pdf"))
    return plt
end

function process_princetonlib_data()
    data = DataFrames.DataFrame(
        [JSON.parse(line) for line in readlines("data/gms.jsonl")],
    )
    data[data[!, :cost] .=== nothing, :cost] .= NaN
    success_df = DataFrames.combine(
        DataFrames.groupby(data, [:file]),
        :cost => (x -> maximum(x) - minimum(x) < 1e-6) => :is_success,

    )
    data = DataFrames.leftjoin(data, success_df, on = :file)
    data = DataFrames.combine(
        DataFrames.groupby(data, [:file, :backend]),
        :time_callbacks => Statistics.mean,
        :time_solve => Statistics.mean,
        :cost => Statistics.mean,
    )
    return data
end

function plot_princetonlib(data, field)
    labels = ["ASL", "CasADi", "Gravity", "JuMP", "SymbolicAD"]
    styles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    colors = Plots.cgrad(:darktest, 5, categorical = true)
    plt = Plots.plot(
        xaxis = :log10,
        legend = :bottomright,
        xlabel = "Time (s)",
        ylabel = "Proportion of problems solved",
        size = (450, 300),
        ylims = (0, 1),
    )
    for label in labels
        df = filter(row -> row.backend == label, data)
        if size(df, 1) == 0
            continue
        end
        full_rows = size(df, 1)
        filter!(row -> !isnan(row.cost_mean), df)
        x = sort(df[!, field])
        if maximum(x) < 1e-6
            continue
        end
        y = (1:size(df, 1)) ./ full_rows
        label = df[1, :backend]
        i = findfirst(isequal(label), labels)
        Plots.plot!(
            plt,
            x,
            y;
            label = label,
            linewidth = 2,
            color = colors[i],
            linestyle = styles[i],
        )
    end
    return plt
end

plot_comparisons(
    "opf-solve.md",
    [
        "JuMP NL" => "ASL",
        "CasADi" => "CasADi",
        "OPOMO" => "Gravity",
        "JuMP" => "JuMP",
        "JuMP SymAD" => "SymbolicAD",
    ],
)


plot_comparisons(
    "opf-callbacks.md",
    [
        "JuMP" => "JuMP",
        "JuMP SymAD" => "SymbolicAD",
    ],
)

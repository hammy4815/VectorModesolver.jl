using CairoMakie

export plot_mode_fields

function _normalize_fields(fields::AbstractVector)
    max_val = maximum(maximum.(fields))
    min_val = minimum(minimum.(fields))
    abs_max = max(max_val,abs(min_val))
    vmax = abs_max
    vmin = -abs_max
    return vmin, vmax
end

"""
    plot_mode_fields(mode_data::Mode; kwargs...)

Plots the six-component, mode-field profiles.
"""
function plot_mode_fields(mode_data::Mode; normalize_fields::Bool=true, kwargs...)
    f = Figure()
    
    labels = [L"E_x", L"Ey", L"Ez"]
    fields = real.([mode_data.Ex, mode_data.Ey, mode_data.Ez])
    vmin, vmax = _normalize_fields(fields)
    for k in 1:3
        ax = Axis(f[2, k],
            title = labels[k],
            xlabel = "X",
            ylabel = "Y",
            aspect=DataAspect(),
        )
        hm = heatmap!(ax, mode_data.x, mode_data.y, fields[k], colormap=:bluesreds, colorrange = (vmin, vmax))
        if k == 1
            Colorbar(f[1, 1:3], hm; vertical=false, label=L"\Re(\textbf{E})")
        end
    end

    labels = [L"Hx", L"Hy", L"Hz"]
    fields = real.([mode_data.Hx, mode_data.Hy, mode_data.Hz])
    vmin, vmax = _normalize_fields(fields)
    for k in 1:3
        ax = Axis(f[4, k],
            title = labels[k],
            xlabel = "X",
            ylabel = "Y",
            aspect=DataAspect(),
        )
        hm = heatmap!(ax, mode_data.x, mode_data.y, fields[k], colormap=:bluesreds, colorrange = (vmin, vmax))
        if k == 1
            Colorbar(f[3, :], hm; vertical=false, label=L"\Re(\textbf{H})")
        end
    end

    resize_to_layout!(f)

    return f
end
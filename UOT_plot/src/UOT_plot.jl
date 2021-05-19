module UOT_plot

using LinearAlgebra, Statistics, PyPlot
using UOT, UOT_toolbox

function plot_points_sources(bounds::Array{Array{Float64,1},1}, x::Array{Array{Float64,1},1}, y::Array{Array{Float64,1},1}, a::Array{Float64, 1},
                                            b::Array{Float64, 1})
    n, m = length(x), length(y);
    figure(figsize = (5, 5))
    for k in 1:n
        ms = 10.0*a[k]/(mean(a)+mean(b));
        if k < n
            plot([x[k][1]], [x[k][2]], linestyle = "none", marker = ".", ms = ms, color = "red")
        else
            plot([x[k][1]], [x[k][2]], linestyle = "none", marker = ".", ms = ms, color = "red", label = L"$\mu_{GT}$")
        end
    end
    for k in 1:m
        ms = 10.0*b[k]/(mean(a)+mean(b));
        if k < m
            plot([y[k][1]], [y[k][2]], linestyle = "none", marker = ".", ms = ms, color = "green")
        else
            plot([y[k][1]], [y[k][2]], linestyle = "none", marker = ".", ms = ms, color = "green", label = L"$\mu_{est}$")
        end
    end

    ax = gca();
    ax.set_xlim([bounds[1][1], bounds[2][1]]);
    ax.set_ylim([bounds[1][2], bounds[2][2]]);
    title(L"$\mu$ and $\nu_{\varepsilon, j}$ with $ε=0.2$")
    legend()
    tight_layout()
    savefig("figures/test.pdf", dpi = 100)
    show()
end

function plot_points_sources_uot(x::Array{Array{Float64,1},1}, y::Array{Array{Float64,1},1}, a::Array{Float64,1},
                                            b::Array{Float64, 1}, bounds::Array{Array{Float64,1},1}, uot_out::UOT.UOT_output)
    n, m = length(x), length(y);
    figure(figsize = (5, 5))
    for k in 1:n
        ms = 25.0*a[k]/(mean(a)+mean(b));
        if k < n
            plot([x[k][1]], [x[k][2]], linestyle = "none", marker = ".", ms = ms, color = "red")
        else
            plot([x[k][1]], [x[k][2]], linestyle = "none", marker = ".", ms = ms, color = "red", label = L"$\mu_{GT}$")
        end
    end
    for k in 1:m
        ms = 25.0*b[k]/(mean(a)+mean(b));
        if k < m
            plot([y[k][1]], [y[k][2]], linestyle = "none", marker = ".", ms = ms, color = "green")
        else
            plot([y[k][1]], [y[k][2]], linestyle = "none", marker = ".", ms = ms, color = "green", label = L"$\mu_{est}$")
        end
    end
    #
    for i in 1:n
        for j in 1:m
            if uot_out.P[i, j] > 1e-8
                lx, ly = UOT_toolbox.line(x[i], y[j]);
                lw = 2*sum(uot_out.P.>0.0)*uot_out.P[i,j]*uot_out.C[i,j]/uot_out.cost;
                if lw > 0.0
                    plot(lx, ly, linestyle = "-", color = "black", lw = lw)
                else
                    println("### error ###")
                end
            end
        end
    end
    #
    for i in 1:n
        if abs(uot_out.mass_variation[i]) > 1e-8
            ms = 25.0*abs(uot_out.mass_variation[i])/(mean(a)+mean(b));
            if uot_out.mass_variation[i] > 0.0
                plot([x[i][1]], [x[i][2]], linestyle = "none", marker = "+", ms = ms, color = "black")
            else
                plot([x[i][1]], [x[i][2]], linestyle = "none", marker = "x", ms = ms, color = "black")
            end
        end
    end
    for j in 1:m
        if abs(uot_out.mass_variation[n + j]) > 1e-8
            ms = 25.0*abs(uot_out.mass_variation[n + j])/(mean(a)+mean(b));
            if uot_out.mass_variation[n + j] > 0.0
                plot([y[j][1]], [y[j][2]], linestyle = "none", marker = "+", ms = ms, color = "black")
            else
                plot([y[j][1]], [y[j][2]], linestyle = "none", marker = "x", ms = ms, color = "black")
            end
        end
    end

    ax = gca();
    ax.set_xlim([bounds[1][1], bounds[2][1]]);
    ax.set_ylim([bounds[1][2], bounds[2][2]]);
    title(string(L"Flat metric between $\mu_{est}$ to $\mu_{GT}$, cost=", round(uot_out.cost; digits = 5)))
    legend()
    tight_layout()
    # savefig("figures/test.pdf", dpi = 100)
    show()
end

end

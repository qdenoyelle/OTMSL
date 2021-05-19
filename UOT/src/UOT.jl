module UOT

using LinearAlgebra
using JuMP
# solver (GLPK could be considered but slower)
using Mosek, MosekTools

function euclidian(x::Array{Float64,1}, y::Array{Float64,1})
    return norm(x - y);
end

mutable struct UOT_output
    cost::Float64
    # matrix of couplings
    P::Array{Float64,2}
    # matrix of costs
    C::Array{Float64,2}
    # vector indicating creation and destruction of mass
    mass_variation::Array{Float64,1}
end

function init_UOT_output(n::Int64 = 1, m::Int64 = 1)
    return UOT_output(Inf, zeros(n, m), zeros(n, m), zeros(n + m));
end

"""
    compute_flat_metric(x, y, a, b, lambda, coupling_mass_variation = true, ground_metric = euclidian)
Compute the flat metric between two point sources distributions given by the locations `x`
with weights `a` and the locations `y` with weights `b`, for the threshold `lambda`. If
`coupling_mass_variation` is set to true then it also computes the coupling and mass variation
between the two distributions.
# Example
```jldoctest
julia> x, y, a, b = [[0.0, 0.0]], [[1.0, 1.0]], [1.0], [0.5];
julia> UOT.compute_flat_metric(x, y, a, b, 1.0)
UOT.UOT_output(1.2071067811865477, [0.5], [1.41421], [0.5, -2.49401e-12])
```
"""
function compute_flat_metric(x::Array{Array{Float64,1},1}, y::Array{Array{Float64,1},1}, a::Array{Float64,1}, b::Array{Float64,1}, lambda::Float64, coupling_mass_variation::Bool = true, ground_metric::Function = euclidian)
    n, m = length(x), length(y);
    out = init_UOT_output(n, m);
    c = cat(a, -b, dims = 1);
    p = n*m;

    # model = Model(with_optimizer(GLPK.Optimizer));
    model = Model(with_optimizer(Mosek.Optimizer, MSK_IPAR_LOG = 0)); # MSK_IPAR_NUM_THREADS number of threads https://www.gams.com/latest/docs/S_MOSEK.html

    @variable(model, f[1:n + m]);
    i = 1;
    for k in 1:n
        for l in 1:m
            e = ground_metric(x[k], y[l]);
            constr = @constraint(model, -e <= f[k] - f[l + n] <= e)
            set_name(constr, "lip_constraints_$(i)")
            i += 1;
        end
    end
    @constraint(model, bound_constraints, -lambda .<= f .<= lambda);

    @objective(model, Max, c'*f);
    JuMP.optimize!(model);
    out.cost = JuMP.objective_value(model);
    # println(termination_status(model))
    # println(JuMP.solve_time(model))

    # compute coupling and cost matrices from primal and dual solutions
    if coupling_mass_variation
        dual = zeros(p);
        for i in 1:p
            dual[i] = -JuMP.dual.(constraint_by_name(model, "lip_constraints_$(i)"));
        end
        out.mass_variation = -JuMP.dual.(bound_constraints);
        f_opt_value = JuMP.value.(f);

        k = 1;
        for i in 1:n
            for j in 1:m
                out.C[i, j] = f_opt_value[i] - f_opt_value[j + n];
                # to prevent numerical approximations
                if abs(dual[k]) > 1e-8
                    out.P[i, j] = dual[k];
                end
                k += 1;
            end
        end
    end
    return out;
end

end # module

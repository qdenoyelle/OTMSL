module UOT_toolbox

using LinearAlgebra

function line(x1::Array{Float64, 1}, x2::Array{Float64, 1})
    t = collect(LinRange(0, 1, 100));
    return [(1-tt)*x1[1] + tt*x2[1] for tt in t], [(1-tt)*x1[2] + tt*x2[2] for tt in t];
end

end

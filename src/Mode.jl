@with_kw struct Mode
    λ::Float64
    neff::Float64
    x::Array{Float64}
    y::Array{Float64}
    Ex::Array{ComplexF64}
    Ey::Array{ComplexF64}
    Ez::Array{ComplexF64}
    Hx::Array{ComplexF64}
    Hy::Array{ComplexF64}
    Hz::Array{ComplexF64}
end
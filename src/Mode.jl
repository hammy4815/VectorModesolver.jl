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

Base.:+(m1::Mode, m2::Mode) = Mode(m1.λ, m1.neff, m1.x, m1.y, 
                                    m1.Ex .+ m2.Ex, m1.Ey .+ m2.Ey, m1.Ez .+ m2.Ez,
                                    m1.Hx .+ m2.Hx, m1.Hy .+ m2.Hy, m1.Hz .+ m2.Hz)
Base.:-(m1::Mode, m2::Mode) = Mode(m1.λ, m1.neff, m1.x, m1.y,
                                    m1.Ex .- m2.Ex, m1.Ey .- m2.Ey, m1.Ez .- m2.Ez,
                                    m1.Hx .- m2.Hx, m1.Hy .- m2.Hy, m1.Hz .- m2.Hz)
Base.:*(α::N, m::Mode) where {N <: Number} = Mode(m.λ, m.neff, m.x, m.y,
                                        α .* m.Ex, α .* m.Ey, α .* m.Ez,
                                        α .* m.Hx, α .* m.Hy, α .* m.Hz)
Base.:*(m::Mode, α::N) where {N <: Number} = Mode(m.λ, m.neff, m.x, m.y,
                                        m.Ex .* α, m.Ey .* α, m.Ez .* α,
                                        m.Hx .* α, m.Hy .* α, m.Hz .* α)
Base.:/(m::Mode, α::N) where {N <: Number} = Mode(m.λ, m.neff, m.x, m.y,
                                        m.Ex ./ α, m.Ey ./ α, m.Ez ./ α,
                                        m.Hx ./ α, m.Hy ./ α, m.Hz ./ α)
    

integrate(x::Array{Float64}, y::Array{Float64}, f::Array{ComplexF64}) = trapz((y,x), f)

"""
"""
poyntingvector(E::Tuple{Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}}, 
            H::Tuple{Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}}) = (
            E[2] .* conj.(H[3]) - E[3] .* conj.(H[2]),
            E[3] .* conj.(H[1]) - E[1] .* conj.(H[3]),
            E[1] .* conj.(H[2]) - E[2] .* conj.(H[1]))
poyntingvector(mode::Mode) = poyntingvector((mode.Ex, mode.Ey, mode.Ez), (mode.Hx, mode.Hy, mode.Hz))

"""
"""
innerproduct(x::Array{Float64}, y::Array{Float64}, 
            E::Tuple{Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}}, 
            H::Tuple{Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}}) = 
            integrate(x, y, poyntingvector(E, H)[3])
innerproduct(mode1::Mode, mode2::Mode) = 
            innerproduct(mode1.x, mode1.y, (mode1.Ex, mode1.Ey, mode1.Ez), (mode2.Hx, mode2.Hy, mode2.Hz))
innerproduct(E::Tuple{Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}}, mode::Mode) =
            innerproduct(mode.x, mode.y, E, H)
innerproduct(mode::Mode, H::Tuple{Array{ComplexF64},Array{ComplexF64},Array{ComplexF64}}) =
            innerproduct(mode.x, mode.y, E, H)

"""
"""
power(x::Array{Float64}, y::Array{Float64},
            E::Tuple{Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}}, 
            H::Tuple{Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}}) = 
            real(innerproduct(x, y, E, H))
power(mode::Mode) = real(innerproduct(mode, mode))

"""
"""
normalize(m::Mode) = m / sqrt(power(m))
function normalize!(m::Mode)
    normfactor = sqrt(power(m))
    m.Ex ./= normfactor 
    m.Ey ./= normfactor
    m.Ez ./= normfactor
    m.Hx ./= normfactor
    m.Hy ./= normfactor
    m.Hz ./= normfactor
    return m
end
"""
"""
modecoeffs(E::Tuple{Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}}, 
            H::Tuple{Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}}, 
            mode::Mode) = (
                innerproduct(mode.x, mode.y, E, (mode.Hx, mode.Hy, mode.Hz)) / power(mode),
                innerproduct(mode.x, mode.y, (mode.Ex, mode.Ey, mode.Ez), H) / power(mode)
            )

function powercoupling(E::Tuple{Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}}, 
                        H::Tuple{Array{ComplexF64}, Array{ComplexF64}, Array{ComplexF64}},
                        mode::Mode)
    a, bstar = modecoeffs(E,H,mode)
    modeip = innerproduct(mode, mode)
    fieldip = innerproduct(mode.x, mode.y, E, H)
    return real(a * bstar * modeip) / real(fieldip)
end
powercoupling(mode1::Mode, mode2::Mode) = powercoupling((mode1.Ex, mode1.Ey, mode1.Ez), 
                                                        (mode1.Hx, mode1.Hy, mode1.Hz), mode2)
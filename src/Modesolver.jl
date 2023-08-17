struct VectorialModesolver
    λ::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    ϵ
end


function assemble(ms::VectorialModesolver)
end

function getHz(Hx, Hy, x, y, β)
    # Init field
    Hz = zeros(ComplexF64, (size(Hx,1)-1,size(Hx,2)-1))
    diffx = x[2:end] .- x[1:end-1]
    diffy = y[2:end] .- y[1:end-1]

    # Get field
    for (i, dx) in enumerate(diffx)
        for (j, dy) in enumerate(diffy)
            Hz[i, j] += (
                (Hx[i+1, j+1] + Hx[i+1, j] - Hx[i, j+1] - Hx[i, j]) / dx + # ∂Hx/∂x
                (Hy[i+1, j+1] + Hy[i, j+1] - Hy[i+1, j] - Hy[i, j]) / dy   # ∂Hy/∂y
            ) / (1im * β)
        end
    end

    # Interpolate field
    itp = interpolate((x, y), Hz, Gridded(Linear()))
    itpHz = zeros(ComplexF64, (size(Hx,1), size(Hx,2)))
    for (i, xv) in enumerate(x)
        for (j, yv) in enumerate(y)
            itpHz[i, j] = itp(xv, yv)
        end
    end

    return itpHz
end

function getE(Hx, Hy, Hz, x, y, β, ω, ϵ)

    # Init Fields
    Ex = zeros(ComplexF64, (size(Hx,1)-1, size(Hx,2)-1))
    Ey = zeros(ComplexF64, (size(Hy,1)-1, size(Hy,2)-1))
    Ez = zeros(ComplexF64, (size(Hz,1)-1, size(Hz,2)-1))
    diffx = x[2:end] .- x[1:end-1]
    diffy = y[2:end] .- y[1:end-1]

    # Get Fields
    # Dx = 1 / (jω) * (∂yHz - ∂zHy) = [∂yHz / (jω)] + (β/ω Hy) 
    # Dy = 1 / (jω) * (∂zHx - ∂xHz) = - [∂zHx / (jω)] - (β/ω Hx)
    # Dz = 1 / (jω) * (∂xHy - ∂yHx) = [∂xHy / (jω)] - [∂yHx / (jω)]
    for (i, dx) in enumerate(diffx)
        for (j, dy) in enumerate(diffy)
            # Get ϵ
            ϵxx, ϵxy, ϵyx, ϵyy, ϵzz = ϵ(x[i], y[j])
            detϵ = ϵxx * ϵyy - ϵxy * ϵyx

            # Build D
            Dx = (
                (Hz[i, j+1] + Hz[i+1, j+1] - Hz[i, j] - Hz[i+1, j]) 
                / (2im * ω * dy) +  # ∂Hz/∂y
                (Hy[i, j] + Hy[i, j+1] + Hy[i+1, j] + Hy[i+1, j+1]) 
                * (0.25 * β/ω)        # -∂Hy/∂z
            )
            Dy = (
                (Hx[i, j] + Hx[i, j+1] + Hx[i+1, j] + Hx[i+1, j+1]) 
                * (-0.25 * β/ω)       # -∂Hx/∂z
                - (Hz[i+1, j] + Hz[i+1, j+1] - Hz[i, j+1] - Hz[i, j]) 
                / (2im * ω * dx)    # ∂Hz/∂x
            )
            Dz = (
                (Hy[i+1, j] + Hy[i+1, j+1] - Hy[i, j+1] - Hy[i, j]) 
                / (2im * ω * dx) +  # ∂Hy/∂x
                (Hx[i, j+1] + Hx[i+1, j+1] - Hx[i, j] - Hx[i+1, j]) 
                / (2im * ω * dy)    # ∂Hx/∂y
            )

            # Get E = ϵ⁻¹ D
            Ex[i, j] = (ϵyy * Dx - ϵxy * Dy) / detϵ
            Ey[i, j] = (ϵxx * Dy - ϵyx * Dx) / detϵ
            Ez[i, j] = Dz / ϵzz
        end
    end

    # Interpolate Fields
    itpEx = interpolate((x, y), Ex, Gridded(Linear()))
    itpEy = interpolate((x, y), Ey, Gridded(Linear()))
    itpEz = interpolate((x, y), Ez, Gridded(Linear()))
    retEx = zeros(ComplexF64, (size(Ex,1), size(Ex,2)))
    retEy = zeros(ComplexF64, (size(Ey,1), size(Ey,2)))
    retEz = zeros(ComplexF64, (size(Ez,1), size(Ez,2)))
    for (i, xv) in enumerate(x)
        for (j, yv) in enumerate(y)
            retEx[i, j] = itpEx(xv, yv)
            retEy[i, j] = itpEy(xv, yv)
            retEz[i, j] = itpEz(xv, yv)
        end
    end

    return retEx, retEy, retEz
end

function solve(ms::VectorialModesolver, nev::Int, tol::Float64, ncv=nothing, sigma=nothing)

    # Get matrix
    A = assemble(ms)

    # Solve eigenvalues
    ncv = ifelse(isnothing(ncv), 10*nev, ncv)
    if isnothing(sigma)
      β²s, ϕs = eigs(A, nev=nev, ncv=ncv, tol=tol, which=:LR)
    else
      β²s, ϕs = eigs(A, nev=nev, ncv=ncv, tol=tol, which=:LM, sigma=sigma)
    end

    # Compile Modes
    modes = Vector{Mode}()
    k = ω = 2 * π / ms.λ
    nx, ny = length(ms.x), length(ms.y)
    for (i, (β², ϕ)) in enumerate(zip(β²s, ϕs))

        # Extract Hx, Hy
        Hx = reshape(ϕ[1:nx*ny, i], (nx, ny))
        Hy = reshape(ϕ[nx*ny+1:end, i], (nx, ny))

        # Get Hz, neff, Ex, Ey, Ez
        β = √(β²)
        Hz = getHz(Hx, Hz, x, y, β)
        neff = β / k
        Ex, Ey, Ez = getE(Hx, Hy, Hz, x, y, β, ω, ms.ϵ)

        # Push Field
        push!(modes, Mode(ms.λ, neff, Ex, Ey, Ez, Hx, Hy, Hz))
    end

    # Sort modes
    sort!(modes, by=m->m.neff)

    return modes
end
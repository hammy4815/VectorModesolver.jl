struct εtype <: Function
end

@with_kw struct VectorialModesolver
    λ::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    boundary::Tuple{Int,Int,Int,Int}
    ε::εtype
end

function assemble(ms::VectorialModesolver)

    # Initialize matrix and vectors
    λ = ms.λ
    k = 2π / λ
    x = ms.x
    y = ms.y
    nx = length(x)
    ny = length(y)
    A = spzeros(Float64, 2 * nx * ny, 2 * nx * ny)
    diffx::Vector{Float64} = x[2:end] .- x[1:end-1]
    diffy::Vector{Float64} = y[2:end] .- y[1:end-1]
    diffx = reshape(vcat(diffx[1], diffx, diffx[end]), :)
    diffy = reshape(vcat(diffy[1], diffy, diffy[end]), :)
    xc::Vector{Float64} = (x[1:end-1] + x[2:end]) ./ 2
    yc::Vector{Float64} = (y[1:end-1] + y[2:end]) ./ 2
    xc = reshape(vcat(xc[1], xc, xc[end]), :)
    yc = reshape(vcat(yc[1], yc, yc[end]), :)

    # Loop through all points
    for i ∈ 1:nx
        for j ∈ 1:ny # Tridiagonal

            # Get ε, nesw
            n = diffy[j+1]
            s = diffy[j]
            e = diffx[i+1]
            w = diffx[i]
            εxx1::Float64, εxy1::Float64, εyx1::Float64, εyy1::Float64, εzz1::Float64 = ms.ε(xc[i], yc[j+1])
            εxx2::Float64, εxy2::Float64, εyx2::Float64, εyy2::Float64, εzz2::Float64 = ms.ε(xc[i], yc[j])
            εxx3::Float64, εxy3::Float64, εyx3::Float64, εyy3::Float64, εzz3::Float64 = ms.ε(xc[i+1], yc[j])
            εxx4::Float64, εxy4::Float64, εyx4::Float64, εyy4::Float64, εzz4::Float64 = ms.ε(xc[i+1], yc[j+1])
            ns21 = n * εyy2 + s * εyy1
            ns34 = n * εyy3 + s * εyy4
            ew14 = e * εxx1 + w * εxx4
            ew23 = e * εxx2 + w * εxx3

            # Eq 21
            axxn = (
                (2 * εyy4 * e - εyx4 * n) * (εyy3 / εzz4) / ns34
                + (2 * εyy1 * w + εyx1 * n) * (εyy2 / εzz1) / ns21
            ) / (n * (e + w))
            # Eq 22
            axxs = (
                (2 * εyy3 * e + εyx3 * s) * (εyy4 / εzz3) / ns34
                + (2 * εyy2 * w - εyx2 * s) * (εyy1 / εzz2) / ns21
            ) / (s * (e + w))
            # Eq 23 transformed
            ayye = (2 * n * εxx4 - e * εxy4) * εxx1 / εzz4 / e / ew14 / (n + s) + (
                2 * s * εxx3 + e * εxy3
            ) * εxx2 / εzz3 / e / ew23 / (n + s)
            # Eq 24 transformed
            ayyw = (2 * εxx1 * n + εxy1 * w) * εxx4 / εzz1 / w / ew14 / (n + s) + (
                2 * εxx2 * s - εxy2 * w
            ) * εxx3 / εzz2 / w / ew23 / (n + s)
            # Eq 23
            axxe = (
                2 / (e * (e + w))
                + (εyy4 * εyx3 / εzz3 - εyy3 * εyx4 / εzz4) / (e + w) / ns34
            )
            # Eq 24
            axxw = (
                2 / (w * (e + w))
                + (εyy2 * εyx1 / εzz1 - εyy1 * εyx2 / εzz2) / (e + w) / ns21
            )
            # Eq 21 transformed
            ayyn = (
                2 / (n * (n + s))
                + (εxx4 * εxy1 / εzz1 - εxx1 * εxy4 / εzz4) / (n + s) / ew14
            )
            # Eq 22 transformed
            ayys = (
                2 / (s * (n + s))
                + (εxx2 * εxy3 / εzz3 - εxx3 * εxy2 / εzz2) / (n + s) / ew23
            )

            # Eq 25
            axxne = +εyx4 * εyy3 / εzz4 / (e + w) / ns34
            axxse = -εyx3 * εyy4 / εzz3 / (e + w) / ns34

            # Eq 26
            axxnw = -εyx1 * εyy2 / εzz1 / (e + w) / ns21
            axxsw = +εyx2 * εyy1 / εzz2 / (e + w) / ns21

            # Eq 25 transformed
            ayyne = +εxy4 * εxx1 / εzz4 / (n + s) / ew14
            ayyse = -εxy3 * εxx2 / εzz3 / (n + s) / ew23

            # Eq 26 transformed
            ayynw = -εxy1 * εxx4 / εzz1 / (n + s) / ew14
            ayysw = +εxy2 * εxx3 / εzz2 / (n + s) / ew23

            # Eq 27
            axxp = (
                - axxn - axxs - axxe - axxw - axxne - axxse - axxnw - axxsw
                + k^2 * (n + s)
                * (εyy4 * εyy3 * e / ns34 + εyy1 * εyy2 * w / ns21)
                / (e + w)
            )

            # Eq 27 transformed
            ayyp = (
                - ayyn - ayys - ayye - ayyw - ayyne - ayyse - ayynw - ayysw
                + k^2 * (e + w)
                * (εxx1 * εxx4 * n / ew14 + εxx2 * εxx3 * s / ew23)
                / (n + s)
            )

            # Eq 28
            axyn = (
                εyy3 * εyy4 / εzz4 / ns34
                - εyy2 * εyy1 / εzz1 / ns21
                + s * (εyy2 * εyy4 - εyy1 * εyy3) / ns21 / ns34
            ) / (e + w)

            # Eq 29
            axys = (
                εyy1 * εyy2 / εzz2 / ns21
                - εyy4 * εyy3 / εzz3 / ns34
                + n * (εyy2 * εyy4 - εyy1 * εyy3) / ns21 / ns34
            ) / (e + w)

            # Eq 28 transformed
            ayxe = (
                εxx1 * εxx4 / εzz4 / ew14
                - εxx2 * εxx3 / εzz3 / ew23
                + w * (εxx2 * εxx4 - εxx1 * εxx3) / ew23 / ew14
            ) / (n + s)

            # Eq 29 transformed
            ayxw = (
                εxx3 * εxx2 / εzz2 / ew23
                - εxx4 * εxx1 / εzz1 / ew14
                + e * (εxx4 * εxx2 - εxx1 * εxx3) / ew23 / ew14
            ) / (n + s)

            # Eq 30
            axye = (εyy4 * (1 + εyy3 / εzz4) - εyy3 * (1 + εyy4 / εzz4)) / ns34 / (
                e + w
            ) - (
                2 * εyx1 * εyy2 / εzz1 * n * w / ns21
                + 2 * εyx2 * εyy1 / εzz2 * s * w / ns21
                + 2 * εyx4 * εyy3 / εzz4 * n * e / ns34
                + 2 * εyx3 * εyy4 / εzz3 * s * e / ns34
                + 2 * εyy1 * εyy2 * (1.0 / εzz1 - 1.0 / εzz2) * w^2 / ns21
            ) / e / (
                e + w
            ) ^ 2

            # Eq 31
            axyw = (εyy2 * (1 + εyy1 / εzz2) - εyy1 * (1 + εyy2 / εzz2)) / ns21 / (
                e + w
            ) - (
                2 * εyx1 * εyy2 / εzz1 * n * e / ns21
                + 2 * εyx2 * εyy1 / εzz2 * s * e / ns21
                + 2 * εyx4 * εyy3 / εzz4 * n * w / ns34
                + 2 * εyx3 * εyy4 / εzz3 * s * w / ns34
                + 2 * εyy3 * εyy4 * (1.0 / εzz3 - 1.0 / εzz4) * e^2 / ns34
            ) / w / (
                e + w
            ) ^ 2

            # Eq 30 transformed
            ayxn = (εxx4 * (1 + εxx1 / εzz4) - εxx1 * (1 + εxx4 / εzz4)) / ew14 / (
                n + s
            ) - (
                2 * εxy3 * εxx2 / εzz3 * e * s / ew23
                + 2 * εxy2 * εxx3 / εzz2 * w * n / ew23
                + 2 * εxy4 * εxx1 / εzz4 * e * s / ew14
                + 2 * εxy1 * εxx4 / εzz1 * w * n / ew14
                + 2 * εxx3 * εxx2 * (1.0 / εzz3 - 1.0 / εzz2) * s^2 / ew23
            ) / n / (
                n + s
            ) ^ 2

            # Eq 31 transformed
            ayxs = (εxx2 * (1 + εxx3 / εzz2) - εxx3 * (1 + εxx2 / εzz2)) / ew23 / (
                n + s
            ) - (
                2 * εxy3 * εxx2 / εzz3 * e * n / ew23
                + 2 * εxy2 * εxx3 / εzz2 * w * n / ew23
                + 2 * εxy4 * εxx1 / εzz4 * e * s / ew14
                + 2 * εxy1 * εxx4 / εzz1 * w * s / ew14
                + 2 * εxx1 * εxx4 * (1.0 / εzz1 - 1.0 / εzz4) * n^2 / ew14
            ) / s / (
                n + s
            ) ^ 2

            # Eq 32
            axyne = +εyy3 * (1 - εyy4 / εzz4) / (e + w) / ns34
            axyse = -εyy4 * (1 - εyy3 / εzz3) / (e + w) / ns34

            # Eq 33
            axynw = -εyy2 * (1 - εyy1 / εzz1) / (e + w) / ns21
            axysw = +εyy1 * (1 - εyy2 / εzz2) / (e + w) / ns21

            # Eq 32 transformed
            ayxne = +εxx1 * (1 - εxx4 / εzz4) / (n + s) / ew14
            ayxse = -εxx2 * (1 - εxx3 / εzz3) / (n + s) / ew23

            # Eq 33 transformed
            ayxnw = -εxx4 * (1 - εxx1 / εzz1) / (n + s) / ew14
            ayxsw = +εxx3 * (1 - εxx2 / εzz2) / (n + s) / ew23

            # Eq 34
            axyp = -(axyn + axys + axye + axyw + axyne + axyse + axynw + axysw) - k^2 * (
                w * (n * εyx1 * εyy2 + s * εyx2 * εyy1) / ns21
                + e * (s * εyx3 * εyy4 + n * εyx4 * εyy3) / ns34
            ) / (e + w)

            # Eq 34 transformed
            ayxp = -(ayxn + ayxs + ayxe + ayxw + ayxne + ayxse + ayxnw + ayxsw) - k^2 * (
                n * (w * εxy1 * εxx4 + e * εxy4 * εxx1) / ew14
                + s * (w * εxy2 * εxx3 + e * εxy3 * εxx2) / ew23
            ) / (n + s)

            # North boundary
            if j == ny
                axxs += ms.boundary[1] * axxn
                axxse += ms.boundary[1] * axxne
                axxsw += ms.boundary[1] * axxnw
                ayxs += ms.boundary[1] * ayxn
                ayxse += ms.boundary[1] * ayxne
                ayxsw += ms.boundary[1] * ayxnw
                ayys -= ms.boundary[1] * ayyn
                ayyse -= ms.boundary[1] * ayyne
                ayysw -= ms.boundary[1] * ayynw
                axys -= ms.boundary[1] * axyn
                axyse -= ms.boundary[1] * axyne
                axysw -= ms.boundary[1] * axynw
            end

            # South boundary
            if j == 1
                axxn += ms.boundary[2] * axxs
                axxne += ms.boundary[2] * axxse
                axxnw += ms.boundary[2] * axxsw
                ayxn += ms.boundary[2] * ayxs
                ayxne += ms.boundary[2] * ayxse
                ayxnw += ms.boundary[2] * ayxsw
                ayyn -= ms.boundary[2] * ayys
                ayyne -= ms.boundary[2] * ayyse
                ayynw -= ms.boundary[2] * ayysw
                axyn -= ms.boundary[2] * axys
                axyne -= ms.boundary[2] * axyse
                axynw -= ms.boundary[2] * axysw
            end

            # East boundary
            if i == nx
                axxw += ms.boundary[3] * axxe
                axxnw += ms.boundary[3] * axxne
                axxsw += ms.boundary[3] * axxse
                ayxw += ms.boundary[3] * ayxe
                ayxnw += ms.boundary[3] * ayxne
                ayxsw += ms.boundary[3] * ayxse
                ayyw -= ms.boundary[3] * ayye
                ayynw -= ms.boundary[3] * ayyne
                ayysw -= ms.boundary[3] * ayyse
                axyw -= ms.boundary[3] * axye
                axynw -= ms.boundary[3] * axyne
                axysw -= ms.boundary[3] * axyse
            end

            # West boundary
            if i == 1
                axxe += ms.boundary[4] * axxw
                axxne += ms.boundary[4] * axxnw
                axxse += ms.boundary[4] * axxsw
                ayxe += ms.boundary[4] * ayxw
                ayxne += ms.boundary[4] * ayxnw
                ayxse += ms.boundary[4] * ayxsw
                ayye -= ms.boundary[4] * ayyw
                ayyne -= ms.boundary[4] * ayynw
                ayyse -= ms.boundary[4] * ayysw
                axye -= ms.boundary[4] * axyw
                axyne -= ms.boundary[4] * axynw
                axyse -= ms.boundary[4] * axysw
            end

            # Construct tensor
            nn = nx * ny

            # Diagonal
            ix = (i - 1) * ny + j # ii
            iy = (i - 1) * ny + j # ii
            A[ix, iy]       = axxp
            A[ix, iy+nn]    = axyp
            A[ix+nn, iy]    = ayxp
            A[ix+nn, iy+nn] = ayyp

            # North
            if (j > 1)
                ix = (i - 1) * ny + j       # n
                iy = (i - 1) * ny + j - 1   # s
                A[ix, iy]       = axxs
                A[ix, iy+nn]    = axys
                A[ix+nn, iy]    = ayxs
                A[ix+nn, iy+nn] = ayys
            end

            # South
            if (j < ny)
                ix = (i - 1) * ny + j     # s
                iy = (i - 1) * ny + j + 1 # n
                A[ix, iy]       = axxn
                A[ix, iy+nn]    = axyn
                A[ix+nn, iy]    = ayxn
                A[ix+nn, iy+nn] = ayyn
            end

            # East
            if (i > 1)
                ix = (i - 1) * ny + j       # e
                iy = (i - 1 - 1) * ny + j   # w
                A[ix, iy]       = axxw
                A[ix, iy+nn]    = axyw
                A[ix+nn, iy]    = ayxw
                A[ix+nn, iy+nn] = ayyw
            end

            # West
            if (i < nx)
                ix = (i - 1) * ny + j         # w
                iy = (i - 1 + 1) * ny + j     # e
                A[ix, iy]       = axxe
                A[ix, iy+nn]    = axye
                A[ix+nn, iy]    = ayxe
                A[ix+nn, iy+nn] = ayye
            end

            # North-East
            if (i > 1 && j > 1)
                ix = (i - 1) * ny + j          # ne
                iy = (i - 1 - 1) * ny + j - 1  # sw
                A[ix, iy]       = axxsw
                A[ix, iy+nn]    = axysw
                A[ix+nn, iy]    = ayxsw
                A[ix+nn, iy+nn] = ayysw
            end

            # South-East
            if (j < ny && i > 1)
                ix = (i - 1) * ny + j          # se
                iy = (i - 1 - 1) * ny + j + 1  # nw
                A[ix, iy]       = axxnw
                A[ix, iy+nn]    = axynw
                A[ix+nn, iy]    = ayxnw
                A[ix+nn, iy+nn] = ayynw
            end

            # South-West
            if (j < ny && i < nx)
                ix = (i - 1) * ny + j          # sw
                iy = (i - 1 + 1) * ny + j + 1  # ne
                A[ix, iy]       = axxne
                A[ix, iy+nn]    = axyne
                A[ix+nn, iy]    = ayxne
                A[ix+nn, iy+nn] = ayyne
            end

            # North-West
            if (j > 1 && i < nx)
                ix = (i - 1) * ny + j          # nw
                iy = (i - 1 + 1) * ny + j - 1  # se
                A[ix, iy]       = axxse
                A[ix, iy+nn]    = axyse
                A[ix+nn, iy]    = ayxse
                A[ix+nn, iy+nn] = ayyse
            end
        end
    end

    return A
end

function getHz(Hx, Hy, x, y, β)
    # Init field
    Hz = zeros(ComplexF64, (size(Hx,1)-1,size(Hx,2)-1))
    diffx = x[2:end] .- x[1:end-1]
    diffy = y[2:end] .- y[1:end-1]

    # Get field
    for (j, dx) in enumerate(diffx)
        for (i, dy) in enumerate(diffy)
            Hz[i, j] = (
                (Hx[i+1, j+1] + Hx[i, j+1] - Hx[i+1, j] - Hx[i, j]) / (2 * dx) + # ∂Hx/∂x
                (Hy[i+1, j+1] + Hy[i+1, j] - Hy[i, j+1] - Hy[i, j]) / (2 * dy)   # ∂Hy/∂y
            ) / (1im * β)
        end
    end

    # Interpolate field
    xc = 0.5 .* (x[2:end] .+ x[1:end-1])
    yc = 0.5 .* (y[2:end] .+ y[1:end-1])
    itp = interpolate((yc, xc), Hz, Gridded(Linear()))
    etpf = extrapolate(itp, Flat())
    itpHz = zeros(ComplexF64, (size(Hx,1), size(Hx,2)))
    for (j, xv) in enumerate(x)
        for (i, yv) in enumerate(y)
            itpHz[i, j] = etpf(yv, xv)
        end
    end

    return itpHz
end

function getE(Hx, Hy, Hz, x, y, β, ω, ε)

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
    for (j, dx) in enumerate(diffx)
        for (i, dy) in enumerate(diffy)
            # Get ε
            xc, yc = (x[j] + x[j+1]) / 2, (y[i] + y[i+1]) / 2
            εxx, εxy, εyx, εyy, εzz = ε(xc, yc)
            detε = εxx * εyy - εxy * εyx

            # Build D
            Dx = (
                (Hz[i+1, j] + Hz[i+1, j+1] - Hz[i, j] - Hz[i, j+1]) 
                / (2im * ω * dy) +  # ∂Hz/∂y
                (Hy[i, j] + Hy[i+1, j] + Hy[i, j+1] + Hy[i+1, j+1]) 
                * (0.25 * β/ω)        # -∂Hy/∂z
            )
            Dy = (
                - (Hx[i, j] + Hx[i, j+1] + Hx[i+1, j] + Hx[i+1, j+1]) 
                * (0.25 * β/ω)       # -∂Hx/∂z
                - (Hz[i+1, j+1] + Hz[i, j+1] - Hz[i+1, j] - Hz[i, j]) 
                / (2im * ω * dx)    # ∂Hz/∂x
            )
            Dz = (
                (Hy[i, j+1] + Hy[i+1, j+1] - Hy[i+1, j] - Hy[i, j]) 
                / (2im * ω * dx) +  # ∂Hy/∂x
                (Hx[i+1, j] + Hx[i+1, j+1] - Hx[i, j] - Hx[i, j+1]) 
                / (-2im * ω * dy)    # ∂Hx/∂y
            )

            # Get E = ε⁻¹ D
            Ex[i, j] = (εyy * Dx - εxy * Dy) / detε
            Ey[i, j] = (εxx * Dy - εyx * Dx) / detε
            Ez[i, j] = Dz / εzz
        end
    end

    # Interpolate Fields
    xc = 0.5 .* (x[2:end] .+ x[1:end-1])
    yc = 0.5 .* (y[2:end] .+ y[1:end-1])
    itpEx = interpolate((yc, xc), Ex, Gridded(Linear()))
    etpEx = extrapolate(itpEx, Flat())
    itpEy = interpolate((yc, xc), Ey, Gridded(Linear()))
    etpEy = extrapolate(itpEy, Flat())
    itpEz = interpolate((yc, xc), Ez, Gridded(Linear()))
    etpEz = extrapolate(itpEz, Flat())
    retEx = zeros(ComplexF64, (size(Ex,1)+1, size(Ex,2)+1))
    retEy = zeros(ComplexF64, (size(Ey,1)+1, size(Ey,2)+1))
    retEz = zeros(ComplexF64, (size(Ez,1)+1, size(Ez,2)+1))
    for (j, xv) in enumerate(x)
        for (i, yv) in enumerate(y)
            retEx[i, j] = etpEx(yv, xv)
            retEy[i, j] = etpEy(yv, xv)
            retEz[i, j] = etpEz(yv, xv)
        end
    end

    return retEx, retEy, retEz
end

function solve(A::SparseMatrixCSC, ms::VectorialModesolver, nev::Int, tol::Float64, ncv=nothing, sigma=nothing, normalize=true)

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
    for (i, β²) in enumerate(β²s)

        # Extract Hx, Hy
        Hx = reshape(ϕs[1:nx*ny, i], (ny, nx))
        Hy = reshape(ϕs[nx*ny+1:end, i], (ny, nx))

        # Get Hz, neff, Ex, Ey, Ez
        β = √(β²)
        Hz = getHz(Hx, Hy, ms.x, ms.y, β)
        neff = β / k
        Ex, Ey, Ez = getE(Hx, Hy, Hz, ms.x, ms.y, β, ω, ms.ε)

        # Push Field
        push!(modes, 
            Mode( 
                λ = ms.λ, 
                neff = neff, 
                x = ms.x, 
                y = ms.y, 
                Ex = Ex, 
                Ey = Ey, 
                Ez = Ez, 
                Hx = Hx, 
                Hy = Hy, 
                Hz = Hz,
                )
            )
    end

    # Sort modes
    sort!(modes, by=m->-m.neff)
    if normalize
        for mode in modes
            normalize!(mode)
        end
    end

    return modes
end

solve(ms::VectorialModesolver, nev::Int, tol::Float64, ncv=nothing, sigma=nothing) = 
    solve(assemble(ms), ms, nev, tol, ncv, sigma)
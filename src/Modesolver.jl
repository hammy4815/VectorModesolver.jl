@with_kw struct Material{T<Number}
    ε_xx::T
    ε_xy::T
    ε_yx::T
    ε_yy::T
    ε_zz::T
end

@with_kw struct Geometry
    shapes::AbstractVector{<:Shape}
    materials::AbstractVector{Material}
end

@with_kw struct VectorialModesolver
    λ::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    boundary::Tuple{Int,Int,Int,Int}
    ϵ::Union{Function, AbstractVector{Geometry}}
end

"""
    get_ε(ε::Function, point::AbstractVector{<:Real})

TBW
"""
get_ε(ε::Function, x::Real, y::Real) = ε(x,y)

"""
    get_ε(ε::AbstractVector{<:Shape}, point::AbstractVector{<:Real})

TBW
"""
function get_ε(ε::AbstractVector{Geometry}, x::Real, y::Real)
    idx = findfirst((x,y), ε.shapes)
    mat = ε.shapes[idx]
    return (mat.ε_xx, mat.ε_xy, mat.ε_yx, mat.ε_yy, mat.ε_zz)
end

"""
    assemble(ms::VectorialModesolver)

TBW
"""
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

            # Get ϵ, nesw
            n = diffy[j+1]
            s = diffy[j]
            e = diffx[i+1]
            w = diffx[i]
            
            ϵxx1::Float64, ϵxy1::Float64, ϵyx1::Float64, ϵyy1::Float64, ϵzz1::Float64 = get_ε(ms.ϵ, xc[i], yc[j+1])
            ϵxx2::Float64, ϵxy2::Float64, ϵyx2::Float64, ϵyy2::Float64, ϵzz2::Float64 = get_ε(ms.ϵ, xc[i], yc[j])
            ϵxx3::Float64, ϵxy3::Float64, ϵyx3::Float64, ϵyy3::Float64, ϵzz3::Float64 = get_ε(ms.ϵ, xc[i+1], yc[j])
            ϵxx4::Float64, ϵxy4::Float64, ϵyx4::Float64, ϵyy4::Float64, ϵzz4::Float64 = get_ε(ms.ϵ, xc[i+1], yc[j+1])

            ns21 = n * ϵyy2 + s * ϵyy1
            ns34 = n * ϵyy3 + s * ϵyy4
            ew14 = e * ϵxx1 + w * ϵxx4
            ew23 = e * ϵxx2 + w * ϵxx3

            # Eq 21
            axxn = (
                (2 * ϵyy4 * e - ϵyx4 * n) * (ϵyy3 / ϵzz4) / ns34
                + (2 * ϵyy1 * w + ϵyx1 * n) * (ϵyy2 / ϵzz1) / ns21
            ) / (n * (e + w))
            # Eq 22
            axxs = (
                (2 * ϵyy3 * e + ϵyx3 * s) * (ϵyy4 / ϵzz3) / ns34
                + (2 * ϵyy2 * w - ϵyx2 * s) * (ϵyy1 / ϵzz2) / ns21
            ) / (s * (e + w))
            # Eq 23 transformed
            ayye = (2 * n * ϵxx4 - e * ϵxy4) * ϵxx1 / ϵzz4 / e / ew14 / (n + s) + (
                2 * s * ϵxx3 + e * ϵxy3
            ) * ϵxx2 / ϵzz3 / e / ew23 / (n + s)
            # Eq 24 transformed
            ayyw = (2 * ϵxx1 * n + ϵxy1 * w) * ϵxx4 / ϵzz1 / w / ew14 / (n + s) + (
                2 * ϵxx2 * s - ϵxy2 * w
            ) * ϵxx3 / ϵzz2 / w / ew23 / (n + s)
            # Eq 23
            axxe = (
                2 / (e * (e + w))
                + (ϵyy4 * ϵyx3 / ϵzz3 - ϵyy3 * ϵyx4 / ϵzz4) / (e + w) / ns34
            )
            # Eq 24
            axxw = (
                2 / (w * (e + w))
                + (ϵyy2 * ϵyx1 / ϵzz1 - ϵyy1 * ϵyx2 / ϵzz2) / (e + w) / ns21
            )
            # Eq 21 transformed
            ayyn = (
                2 / (n * (n + s))
                + (ϵxx4 * ϵxy1 / ϵzz1 - ϵxx1 * ϵxy4 / ϵzz4) / (n + s) / ew14
            )
            # Eq 22 transformed
            ayys = (
                2 / (s * (n + s))
                + (ϵxx2 * ϵxy3 / ϵzz3 - ϵxx3 * ϵxy2 / ϵzz2) / (n + s) / ew23
            )

            # Eq 25
            axxne = +ϵyx4 * ϵyy3 / ϵzz4 / (e + w) / ns34
            axxse = -ϵyx3 * ϵyy4 / ϵzz3 / (e + w) / ns34

            # Eq 26
            axxnw = -ϵyx1 * ϵyy2 / ϵzz1 / (e + w) / ns21
            axxsw = +ϵyx2 * ϵyy1 / ϵzz2 / (e + w) / ns21

            # Eq 25 transformed
            ayyne = +ϵxy4 * ϵxx1 / ϵzz4 / (n + s) / ew14
            ayyse = -ϵxy3 * ϵxx2 / ϵzz3 / (n + s) / ew23

            # Eq 26 transformed
            ayynw = -ϵxy1 * ϵxx4 / ϵzz1 / (n + s) / ew14
            ayysw = +ϵxy2 * ϵxx3 / ϵzz2 / (n + s) / ew23

            # Eq 27
            axxp = (
                - axxn - axxs - axxe - axxw - axxne - axxse - axxnw - axxsw
                + k^2 * (n + s)
                * (ϵyy4 * ϵyy3 * e / ns34 + ϵyy1 * ϵyy2 * w / ns21)
                / (e + w)
            )

            # Eq 27 transformed
            ayyp = (
                - ayyn - ayys - ayye - ayyw - ayyne - ayyse - ayynw - ayysw
                + k^2 * (e + w)
                * (ϵxx1 * ϵxx4 * n / ew14 + ϵxx2 * ϵxx3 * s / ew23)
                / (n + s)
            )

            # Eq 28
            axyn = (
                ϵyy3 * ϵyy4 / ϵzz4 / ns34
                - ϵyy2 * ϵyy1 / ϵzz1 / ns21
                + s * (ϵyy2 * ϵyy4 - ϵyy1 * ϵyy3) / ns21 / ns34
            ) / (e + w)

            # Eq 29
            axys = (
                ϵyy1 * ϵyy2 / ϵzz2 / ns21
                - ϵyy4 * ϵyy3 / ϵzz3 / ns34
                + n * (ϵyy2 * ϵyy4 - ϵyy1 * ϵyy3) / ns21 / ns34
            ) / (e + w)

            # Eq 28 transformed
            ayxe = (
                ϵxx1 * ϵxx4 / ϵzz4 / ew14
                - ϵxx2 * ϵxx3 / ϵzz3 / ew23
                + w * (ϵxx2 * ϵxx4 - ϵxx1 * ϵxx3) / ew23 / ew14
            ) / (n + s)

            # Eq 29 transformed
            ayxw = (
                ϵxx3 * ϵxx2 / ϵzz2 / ew23
                - ϵxx4 * ϵxx1 / ϵzz1 / ew14
                + e * (ϵxx4 * ϵxx2 - ϵxx1 * ϵxx3) / ew23 / ew14
            ) / (n + s)

            # Eq 30
            axye = (ϵyy4 * (1 + ϵyy3 / ϵzz4) - ϵyy3 * (1 + ϵyy4 / ϵzz4)) / ns34 / (
                e + w
            ) - (
                2 * ϵyx1 * ϵyy2 / ϵzz1 * n * w / ns21
                + 2 * ϵyx2 * ϵyy1 / ϵzz2 * s * w / ns21
                + 2 * ϵyx4 * ϵyy3 / ϵzz4 * n * e / ns34
                + 2 * ϵyx3 * ϵyy4 / ϵzz3 * s * e / ns34
                + 2 * ϵyy1 * ϵyy2 * (1.0 / ϵzz1 - 1.0 / ϵzz2) * w^2 / ns21
            ) / e / (
                e + w
            ) ^ 2

            # Eq 31
            axyw = (ϵyy2 * (1 + ϵyy1 / ϵzz2) - ϵyy1 * (1 + ϵyy2 / ϵzz2)) / ns21 / (
                e + w
            ) - (
                2 * ϵyx1 * ϵyy2 / ϵzz1 * n * e / ns21
                + 2 * ϵyx2 * ϵyy1 / ϵzz2 * s * e / ns21
                + 2 * ϵyx4 * ϵyy3 / ϵzz4 * n * w / ns34
                + 2 * ϵyx3 * ϵyy4 / ϵzz3 * s * w / ns34
                + 2 * ϵyy3 * ϵyy4 * (1.0 / ϵzz3 - 1.0 / ϵzz4) * e^2 / ns34
            ) / w / (
                e + w
            ) ^ 2

            # Eq 30 transformed
            ayxn = (ϵxx4 * (1 + ϵxx1 / ϵzz4) - ϵxx1 * (1 + ϵxx4 / ϵzz4)) / ew14 / (
                n + s
            ) - (
                2 * ϵxy3 * ϵxx2 / ϵzz3 * e * s / ew23
                + 2 * ϵxy2 * ϵxx3 / ϵzz2 * w * n / ew23
                + 2 * ϵxy4 * ϵxx1 / ϵzz4 * e * s / ew14
                + 2 * ϵxy1 * ϵxx4 / ϵzz1 * w * n / ew14
                + 2 * ϵxx3 * ϵxx2 * (1.0 / ϵzz3 - 1.0 / ϵzz2) * s^2 / ew23
            ) / n / (
                n + s
            ) ^ 2

            # Eq 31 transformed
            ayxs = (ϵxx2 * (1 + ϵxx3 / ϵzz2) - ϵxx3 * (1 + ϵxx2 / ϵzz2)) / ew23 / (
                n + s
            ) - (
                2 * ϵxy3 * ϵxx2 / ϵzz3 * e * n / ew23
                + 2 * ϵxy2 * ϵxx3 / ϵzz2 * w * n / ew23
                + 2 * ϵxy4 * ϵxx1 / ϵzz4 * e * s / ew14
                + 2 * ϵxy1 * ϵxx4 / ϵzz1 * w * s / ew14
                + 2 * ϵxx1 * ϵxx4 * (1.0 / ϵzz1 - 1.0 / ϵzz4) * n^2 / ew14
            ) / s / (
                n + s
            ) ^ 2

            # Eq 32
            axyne = +ϵyy3 * (1 - ϵyy4 / ϵzz4) / (e + w) / ns34
            axyse = -ϵyy4 * (1 - ϵyy3 / ϵzz3) / (e + w) / ns34

            # Eq 33
            axynw = -ϵyy2 * (1 - ϵyy1 / ϵzz1) / (e + w) / ns21
            axysw = +ϵyy1 * (1 - ϵyy2 / ϵzz2) / (e + w) / ns21

            # Eq 32 transformed
            ayxne = +ϵxx1 * (1 - ϵxx4 / ϵzz4) / (n + s) / ew14
            ayxse = -ϵxx2 * (1 - ϵxx3 / ϵzz3) / (n + s) / ew23

            # Eq 33 transformed
            ayxnw = -ϵxx4 * (1 - ϵxx1 / ϵzz1) / (n + s) / ew14
            ayxsw = +ϵxx3 * (1 - ϵxx2 / ϵzz2) / (n + s) / ew23

            # Eq 34
            axyp = -(axyn + axys + axye + axyw + axyne + axyse + axynw + axysw) - k^2 * (
                w * (n * ϵyx1 * ϵyy2 + s * ϵyx2 * ϵyy1) / ns21
                + e * (s * ϵyx3 * ϵyy4 + n * ϵyx4 * ϵyy3) / ns34
            ) / (e + w)

            # Eq 34 transformed
            ayxp = -(ayxn + ayxs + ayxe + ayxw + ayxne + ayxse + ayxnw + ayxsw) - k^2 * (
                n * (w * ϵxy1 * ϵxx4 + e * ϵxy4 * ϵxx1) / ew14
                + s * (w * ϵxy2 * ϵxx3 + e * ϵxy3 * ϵxx2) / ew23
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
    for (j, dx) in enumerate(diffx)
        for (i, dy) in enumerate(diffy)
            # Get ϵ
            xc, yc = (x[j] + x[j+1]) / 2, (y[i] + y[i+1]) / 2
            ϵxx, ϵxy, ϵyx, ϵyy, ϵzz = ϵ(xc, yc)
            detϵ = ϵxx * ϵyy - ϵxy * ϵyx

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

            # Get E = ϵ⁻¹ D
            Ex[i, j] = (ϵyy * Dx - ϵxy * Dy) / detϵ
            Ey[i, j] = (ϵxx * Dy - ϵyx * Dx) / detϵ
            Ez[i, j] = Dz / ϵzz
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

function solve(A::SparseMatrixCSC, ms::VectorialModesolver, nev::Int, tol::Float64, ncv=nothing, sigma=nothing)

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
        Ex, Ey, Ez = getE(Hx, Hy, Hz, ms.x, ms.y, β, ω, ms.ϵ)

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

    return modes
end

solve(ms::VectorialModesolver, nev::Int, tol::Float64, ncv=nothing, sigma=nothing) = 
    solve(assemble(ms), ms, nev, tol, ncv, sigma)
using Base.Threads

"""
	gaussianblur(I, order)-> img

Two-dimensional separable gaussian blur.
Low pass image filtering, smooth, circular symetry,
respect edges in various directions.
"""
function gaussianblur(I::Array{Float32,2})
    h, w = size(I)
    #iseven(order) && (1 < order < 13) && error("Order must be even [3,11]")
    c = [0.375f0 0.25f0 0.0625f0] # coefficients for order 5
    J = Array{Float32,2}(undef, h, w)

    for j = 1:w, i = 3:h-2
        @inbounds J[i, j] =
            I[i, j] * c[1] +
            (I[i-1, j] + I[i+1, j]) * c[2] +
            (I[i-2, j] + I[i+2, j]) * c[3]
    end
    for i in [1 2 h - 1 h]
        J[i, :] .= I[i, :]
    end

    K = Array{Float32,2}(undef, h, w)
    for j = 3:w-2, i = 1:h
        @inbounds K[i, j] =
            J[i, j] * c[1] +
            (J[i, j-1] + J[i, j+1]) * c[2] +
            (J[i, j-2] + J[i, j+2]) * c[3]
    end
    for j in [1 2 w - 1 w]
        K[:, j] .= J[:, j]
    end

    return K
end

function gaussianblur1(I::Array{Float32,2})
    h, w = size(I)
    #iseven(order) && (1 < order < 13) && error("Order must be even [3,11]")

    c = [6.0f0 4.0f0]
    a = 0.0625f0
    J = Array{Float32,2}(undef, h, w)
    @inbounds for j = 1:w, i = 3:h-2
        J[i, j] =
            (
                I[i, j] * c[1] +
                (I[i-1, j] + I[i+1, j]) * c[2] +
                (I[i-2, j] + I[i+2, j])
            )a
    end
    for i in [1 2 h - 1 h]
        J[i, :] .= I[i, :]
    end

    K = Array{Float32,2}(undef, h, w)
    @inbounds for j = 3:w-2, i = 1:h
        K[i, j] =
            (
                J[i, j] * c[1] +
                (J[i, j-1] + J[i, j+1]) * c[2] +
                (J[i, j-2] + J[i, j+2])
            )a
    end
    for j in [1 2 w - 1 w]
        K[:, j] .= J[:, j]
    end

    return K
end

function gaussianblur2(I::Array{Float32,2})
    h, w = size(I)
    #iseven(order) && (1 < order < 13) && error("Order must be even [3,11]")
    c = [0.375f0 0.25f0 0.0625f0] # coefficients for order 5

    J = Array{Float32,2}(undef, h, w)
    @inbounds for j = 1:w, i = 3:h-2
        J[i, j] = I[i, j] * c[1]
        for k = 1:2
            J[i, j] += (I[i-k, j] + I[i+i, j]) * c[k+1]
        end
    end
    for i in [1 2 h - 1 h]
        J[i, :] .= I[i, :]
    end

    K = Array{Float32,2}(undef, h, w)
    @inbounds for j = 3:w-2, i = 1:h
        K[i, j] = J[i, j] * c[1]
        for k = 1:2
            K[i, j] = (J[i, j-k] + J[i, j+k]) * c[k+1]
        end
    end
    for j in [1 2 w - 1 w]
        K[:, j] .= J[:, j]
    end

    return K
end

function gaussianblur3(I::Array{Float32,2})
    h, w = size(I)
    n = h * w
    c = [0.375f0 0.25f0 0.0625f0]
    h2 = h << 1
    n1 = n + 1
    J = Array{Float32,2}(undef, h, w)
    for i = 3:n-2
        @inbounds J[i] =
            I[i] * c[1] + (I[i-1] + I[i+1]) * c[2] + (I[i-2] + I[i+2]) * c[3]
    end
    for i in [1 2 n - 1 n]
        J[i] = I[i]
    end

    K = Array{Float32,2}(undef, h, w)
    for i = 1+h2:n-h2
        @inbounds K[i] =
            J[i] * c[1] + (J[i-h] + J[i+h]) * c[2] + (J[i-h2] + J[i+h2]) * c[3]
    end

    for i = 1:h2
        @inbounds K[i] = K[i+h2]
    end
    for i = n+1-2h:n
        @inbounds K[i] = K[i-2h]
    end

    return K
end

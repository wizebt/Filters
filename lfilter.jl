# Bilinear filter 5x5
function lfilter(I::Array{Float32,2}) # 5x5
    h, w = size(I)
    n = h * w
    J = Array{Float32,1}(undef, n)

    for i in [1 2 n - 1 n]
        @inbounds J[i] = I[i]
    end
    for i = 3:n-2
        @inbounds J[i] = (I[i] + I[i-1] + I[i+1] + I[i-2] + I[i+2])/5.0f0
    end
    return J

    K = Array{Float32,2}(undef, h, w)
    h2 = h << 1
    for i = 1+h2:n-h2
        @inbounds K[i] = (J[i-h2] + J[i-h] + J[i] + J[i+h] + J[i+h2])/5.0f0
    end
    for i = 1:h2
        @inbounds K[i] = J[i]
    end
    for i = 1+n-h2:n
        @inbounds K[i] = J[i]
    end
    return K
end

using BenchmarkTools
I = ones(Float32, 5, 10);
J = rand(Float32, 41, 41);

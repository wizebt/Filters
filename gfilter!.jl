"""
    gfilter(I::Array{Float32,2}) -> gaussian filtered 

Fast gaussian 5x5 filtering, using separation and kernel folding.
"""
function gfilter(I::Array{Float32,2}) # 5x5
    h, w = size(I)
    n = h * w
    #a, b, c = 0.375f0, 0.25f0, 0.0625f0 # [1 4 6 4 1] / 16
    a, b, c = 0.40262f0, 0.244201f0, 0.0544887f0

    J = Array{Float32,1}(undef, n)
    for i in [1 2 n - 1 n]
        @inbounds J[i] = I[i]
    end
    for i = 3:n-2
        @inbounds J[i] = I[i]a + (I[i-1] + I[i+1])b + (I[i-2] + I[i+2])c
    end

    K = Array{Float32,2}(undef, h, w)
    h2 = h << 1
    for i = 1+h2:n-h2
        @inbounds K[i] = J[i]a + (J[i-h] + J[i+h])b + (J[i-h2] + J[i+h2])c
    end
    for i = 1:h2
        @inbounds K[i] = J[i]
    end
    for i = 1+n-h2:n
        @inbounds K[i] = J[i]
    end

    return K
end

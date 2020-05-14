using Images, ImageView, BenchmarkTools, Base.Threads

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

function gfilter1(I::Array{Float32,2}) # 3x3
    h, w = size(I)
    n = h * w
    a, b = 0.5f0, 0.25f0

    J = Array{Float32, 1}(undef, n)
    # Replicate start and end
    J[1] = I[1]; J[n] = I[n]
    for i = 2:n-1
        @inbounds J[i] = I[i]a + (I[i-1] + I[i+1])b
    end

    K = Array{Float32,2}(undef, h, w)
    for i = 1+h:n-h
        @inbounds K[i] = J[i]a + (J[i-h] + J[i+h])b
    end

    @inbounds for i = 1:h
        K[i] = K[i+h]
        K[n+1-i] = K[n+1-h]
    end

    return K
end

function gfilterwhile(I)
    h, w = size(I)
    n = h * w
    #a, b, c = 0.40262f0, 0.244201f0, 0.0544887f0
    #a, b, i, n = I[1], I[2], 2, length(I)
    a, b = 0.5f0, 0.25f0
    i = 2
    J = Array{Float32,1}(undef, n)
    # Replicate start and end
    J[1] = I[1]; J[n] = I[n]
    while i < n
        @inbounds J[i] = I[i]a + (I[i-1] + I[i+1])b
        i += 1
    end
    K = Array{Float32,2}(undef, h, w)
    i = 2
    while i < n
        @inbounds K[i] = J[i]a + (J[i-h] + J[i+h])b
        i += 1
    end
    @inbounds for i = 1:h
        K[i] = K[i+h]
        K[n+1-i] = K[n+1-h]
    end
    return K
end

function gfilter!(I)
    a, b, c, i, n = I[1], I[2], 2, length(I)
    while i < n
        i1 = i + 1
        c = I[i1]
        I[i] = 0.5b + 0.25 * (a + c)
        a, b, i = b, c, i1
    end

    a, b, c, i = I[1], I[2], 2
    while i < n
        i1 = i + 1
        c = I[i1]
        I[i] = 0.5b + 0.25 * (a + c)
        a, b, i = b, c, i1
    end
    return transpose(I)
end

function gfilteri(I::Array{UInt32,2})
    g = [6 4 1] # / 16
    h, w = size(I)
    n = h * w
    J = Array{UInt32,2}(undef, h, w)

    for i = 3:n-2
        @inbounds J[i] =
            (6 * I[i] + (I[i-1] + I[i+1]) >> 2 + I[i-2] + I[i+2]) >> 4
    end

    # Replicate start and end
    for i = 1:2
        J[i] = J[i+2]
        J[n-i+1] = J[n-i-1]
    end

    K = Array{UInt32,2}(undef, h, w)
    for i = 1+2h:n-2h
        @inbounds K[i] =
            (6 * J[i] + (J[i-h] + J[i+h]) >> 2 + J[i-2h] + J[i+2h]) >> 4
    end

    for i = 1:2h
        @inbounds K[i] = K[n+1-i] = K[i+2h]
    end

    return K
end

function gfilter1!(J::Array{Float32,2}, I::Array{Float32,2})
    a, b, c = 0.375f0, 0.25f0, 0.0625f0 # [1 4 6 4 1] / 16
    h, w = size(I)
    n = h * w

    K = Array{Float32,2}(undef, h, w)
    for i = 3:h-2, j in [1 2 n - 1 n]
        K[i, j] = I[i, j]
    end
    for i = 3:n-2
        @inbounds J[i] = I[i]a + (I[i-1] + I[i+1])b + (I[i-2] + I[i+2])c
    end
    for i = 1+2h:n-2h
        @inbounds K[i] = J[i]a + (J[i-h] + J[i+h])b + (J[i-2h] + J[i+2h])c
    end
    for i = 1:2h
        @inbounds K[i] = K[n+1-i] = K[i+2h]
    end

    return K
end

function t(I)
    a, b, i, n = I[1], I[2], 2, length(I)
    while i < n
        i1 = i + 1
        c = I[i1]
        I[i] = a + 2b + c
        a, b, i = b, c, i1
    end
    J = transpose(I)
    a, b, i = J[1], J[2], 2
    while i < n
        i1 = i + 1
        c = J[i1]
        J[i] = a + 2b + c
        a, b, i = b, c, i1
    end

    return transpose(J)
end

# This thread version is slower for small array
using Base.Threads
function gfilter1(I) # 5x5
    a, b, c = 0.375f0, 0.25f0, 0.0625f0 # [1 4 6 4 1] / 16
    h, w = size(I)
    n = h * w
    z = 0.0f0
    J = Array{Float32,2}(undef, h, w)
    J[[1 2 end - 1 end]] .= z

    @threads for i = 3:n-2
        @inbounds J[i] = I[i]a + (I[i-1] + I[i+1])b + (I[i-2] + I[i+2])c
    end

    K = Array{Float32,2}(undef, h, w)
    @threads for i = 1+2h:n-2h
        @inbounds K[i] = J[i]a + (J[i-h] + J[i+h])b + (J[i-2h] + J[i+2h])c
    end

    @threads for i = 1:2h
        @inbounds K[i] = K[n+1-i] = K[i+2h]
    end

    return K
end

function t!(J::Array{Float32,2}, I::Array{Float32,2})
    a, b, c = 0.375f0, 0.25f0, 0.0625f0
    @threads for i = 3:length(I)-2
        @inbounds J[i] = I[i]a + (I[i-1] + I[i+1])b + (I[i-2] + I[i+2])c
    end
end

"""
	gaussianblur3(I::Array{Float32,2}) -> gaussian filtered image

Gaussian blur 3x3
"""
function gaussianblur3(I::Array{Float32,2})
    # gaussian g = [1 2 1]/4
    h, w = size(I)
    n = h * w

    A = Array{Float32,1}(undef, n)
    for i = 2:n-1
        @inbounds A[i] = 0.5f0 * I[i] + 0.25f0 * (I[i-1] + I[i+1])
    end
    A[1] = A[end] = 0.0f0
    J = Array{Float32,2}(undef, h, w)
    for i = 1+h:n-h
        @inbounds J[i] = 0.5f0 * A[i] + 0.25f0 * (A[i-h] + A[i+h])
    end

    return J
end
"""
	gaussianblur3(I::Array{Float32,2}) -> gaussian filtered image

Gaussian blur 3x3
"""
function gaussianblur3!(J::Array{Float32,2}, I::Array{Float32,2})
    # gaussian g = [1 2 1]/4
    h, w = size(I)
    n = h * w

    A = Array{Float32,1}(undef, n)
    for i = 2:n-1
        A[i] = 0.5f0 * I[i] + 0.25f0 * (I[i-1] + I[i+1])
    end
    A[1] = A[end] = 0.0f0

    for i = 1+h:n-h
        J[i] = 0.5f0 * A[i] + 0.25f0 * (A[i-h] + A[i+h])
    end

    return J
end

"""
	gaussianblur5(I::Array{Float32,2}) -> gaussian filtered image

Gaussian blur 5x5
"""
function gaussianblur5(I::Array{Float32,2})
    # gaussian g = [1 4 6 4 1]/16
    h, w = size(I)
    n = h * w

    A = Array{Float32,1}(undef, n)
    for i = 3:n-2
        @inbounds A[i] =
            0.375f0 * I[i] +
            0.25f0 * (I[i-1] + I[i+1]) +
            0.0625f0 * (I[i-2] + I[i+2])
    end
    A[[1 2 end - 1 end]] .= 0.0f0

    J = Array{Float32,2}(undef, h, w)
    for i = 1+2h:n-2h
        @inbounds J[i] =
            0.375f0 * A[i] +
            0.25f0 * (A[i-h] + A[i+h]) +
            0.0625f0 * (A[i-2h] + A[i+2h])
    end
    return J
end

function bing1(I)
    g = [0.40262f0 0.244201f0 0.0544887f0]
    h, w = size(I)
    J = Array{Float32,2}(undef, h, w)
    n = h * w
    z = 0.0f0

    for i in [1 2 n - 1 n]
        J[i] = z
    end

    for i = 3:n-2
        @inbounds J[i] =
            I[i] * g[1] + (I[i-1] + I[i+1]) * g[2] + (I[i-2] + I[i+2]) * g[3]
    end

    K = Array{Float32,2}(undef, h, w)
    mx, mi = z, 1 + 2h
    for i = 1+2h:n-2h
        @inbounds K[i] =
            J[i] * g[1] + (J[i-h] + J[i+h]) * g[2] + (J[i-2h] + J[i+2h]) * g[3]
        if K[i] > mx
            mx, mi = K[i], i
        end
    end

    return mx, rem(mi - 1, h) + 1, div(mi - 1, h) + 1
end

using ImageFiltering, ImageView
gk1 = Kernel.gaussian(1.0f0);
gk2 = Kernel.gaussian(2.0f0);
gk3 = Kernel.gaussian(3.0f0);
g = collect(gk2);

h, w = 60, 80
I = zeros(Float32, h, w);
n = rand(Float32, h, w);
I = I + 2n;
g = collect(Kernel.gaussian(1.0f0))
I[11:15, 21:25] .+= 40g; # centroid 15,25
g = collect(gk2)
I[41:49, 51:59] .+= 120g; # centroid 45,55
for i = 1:99 I[rand(1:60), rand(1:80)] += 9.0f0 end
imshow(I, name = "Test image")
imshow(median33(I), name="Median filter")
imshow(gfilter(I), name = "gfilter");
imshow(imfilter(I, gk2), name = "imfilter 2");
imshow(imfilter(I, gk1), name = "imfilter 1");

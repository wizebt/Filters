function L(i::Int)
	N=2r+1
	l=Array{Float32,1}(undef,N)
	l[r+1]=1
	for i=1:r
		l[r+1+i]=l[r+1-i] = sinc(pi*(2k/(N-1)-1))
	end
	return l
end

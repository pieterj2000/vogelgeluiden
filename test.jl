# some imports
using Plots
using FFTW

win_len = 50
win = sin.(π * (1:win_len)/(win_len+1)).^2;
hop = 15

plot(win)

norm = win.^2
for i = 1+hop:hop:win_len
    norm[1:end-i+1] += win[i:end].^2
    norm[i:end] += win[1:end-i+1].^2
end

win_hat = win ./ sqrt.(real(norm))
win = win_hat
plot(win)
     

K = 30
signal_len = win_len + (K-1) * hop
x = randn(signal_len);
plot(x)


# declare a complex array to hold the STFT of x
X = Array{Complex,2}(undef, (win_len, K))

# compute windowed fft's
for (i,j) = zip( 1 : hop : (signal_len - win_len + 1), 1 : K )
    X[:,j] = fft(x[i:i + win_len - 1] .* win) / √(win_len) 
end

# definition of Map1, M1
function Map1(x, signal_length, window_length, hop_size, hop_count)
    X = Array{Complex, 2}(undef, (window_length,K)) 
    for (i,j) = zip( 1 : hop_size : (signal_length - window_length + 1), 1 : hop_count )
        X[:,j] = x[i:i + window_length - 1]
    end
    return X
end
M1 = z -> Map1(z, signal_len, win_len, hop, K)

# definition of Map2, M2
function Map2(X, window)
    return X .* window
end

M2 = Z -> Map2(Z, win)

# definition of M3, as an orthogonal operation
M3 = Z -> fft(Z,[1])/√(win_len)

# we declare a new map by cascading these operations.
M = z -> M3( M2( M1(z) ) )
     
Y = M(x);
error = sum(abs.(Y[:] - X[:]).^2)
println("Error : ", error)

M3T = Z -> ifft(Z,[1]) * √(win_len)

M2T = M2

function Map1Xpose(X, signal_len, win_len, hop, K)
    x = zeros(signal_len)im
    for (i,j) = zip( 1 : hop : (signal_len - win_len + 1), 1 : K )
        x[i:i + win_len - 1] += X[:,j]
    end
    return x
end
M1T = Z -> Map1Xpose(Z, signal_len, win_len, hop, K)


inner_product = (x,y) -> sum( (x .* conj(y))[:] )

z = randn(signal_len)+ randn(signal_len)im
Z = randn(win_len,K) + randn(win_len,K)im
in1 = inner_product( M1(z), Z )
in2 = inner_product( z, M1T(Z) )
println("*** Testing transpose for M1 ***")
println("Difference between the two inner products: ", in1 - in2)


Z2 = randn(win_len,K) + randn(win_len,K)im
in1 = inner_product( M2(Z), Z2 )
in2 = inner_product( Z, M2T(Z2) )
println("*** Testing transpose for M2 ***")
println("Difference between the two inner products: ", in1 - in2)


in1 = inner_product( M3(Z), Z2 )
in2 = inner_product( Z, M3T(Z2) )
println("*** Testing transpose for M3 ***")
println("Difference between the two inner products: ", in1 - in2)


MT = Z -> M1T( M2T( M3T(Z) ) )

MTM = z -> MT( M(z) )


x2 = MTM(x) .|> real
plot([x, x2])
plot((x-x2)[100:400])

ax.legend()
     
P = MTM(ones(signal_len)) .|> real
plot(P)
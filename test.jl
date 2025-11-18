# some imports
using Plots
using FFTW

include("stft.jl")
using .STFT

win_len = 50
win = STFT.hann_window(win_len)
hop = 15

win == sin.(π * (1:win_len)/(win_len+1)).^2
plot(win)

norm = win.^2
for i = 1+hop:hop:win_len
    norm[1:end-i+1] += win[i:end].^2
    norm[i:end] += win[1:end-i+1].^2
end

win_hat = win ./ sqrt.(real(norm))
#win = win_hat
plot(win)
     

K = 30
signal_len = win_len + (K-1) * hop
x = randn(signal_len)
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
     
Y = M(x)
Z = STFT.stft(Complex.(x), Complex.(win), hop)
#Z = STFT.rstft(x, win, hop)
error = sum(abs.(Y[:] - X[:]).^2)
error = sum(abs.(Z[:] - X[:]).^2)
println("Error : ", error)
sum(abs2.(Y-Z))

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




z = randn(signal_len)+ randn(signal_len)im
Z = randn(win_len,K) + randn(win_len,K)im
in1 = inner_product( STFT.maaksegments(z, win_len, hop), Z )
in2 = inner_product( z, STFT.overlapadd(Z, hop) )
in2 = inner_product( z, M1T(Z) )
println("*** Testing transpose for M1 ***")
println("Difference between the two inner products: ", in1 - in2)



Z2 = randn(win_len,K) + randn(win_len,K)im
in1 = inner_product( STFT.applyWindow(Z, win), Z2 )
in2 = inner_product( Z, STFT.applyWindow(Z2, win) )
println("*** Testing transpose for M2 ***")
println("Difference between the two inner products: ", in1 - in2)


in1 = inner_product( M3(Z), Z2 )
in2 = inner_product( Z, M3T(Z2) )
println("*** Testing transpose for M3 ***")
println("Difference between the two inner products: ", in1 - in2)


MT = Z -> M1T( M2T( M3T(Z) ) )

MTM = z -> MT( M(z) )

zhun = MT(Y)
zons = STFT.istft(Y, win, hop)
zhun-zons
sum(abs2.(zhun-zons))

MTMons = z -> STFT.istft(STFT.stft(Complex.(z), win, hop) , win, hop)

x2 = MTM(x) #.|> real
interm = STFT.stft(Complex.(x), win, hop)
x2ons = STFT.istft(interm , win, hop)
x2ons = MTMons(x)

abs2.(x2 - x2ons) |> sum
abs2.(x - x2ons) |> sum

x2 = real.(x2)
x2ons = real.(x2ons)

plot([x, x2ons])

MTM(ones(signal_len))
MTMons(ones(signal_len))

win_hatons = STFT.normalize_window(win, hop)
plot(win)
plot!(win_hatons)
z = randn(signal_len)+ randn(signal_len)im
Z = randn(win_len,K) + randn(win_len,K)im
in1 = inner_product( STFT.stft(z, win, hop), Z )
in2 = inner_product( z, STFT.istft(Z , win, hop) )
println("*** Testing transpose for M3 ***")
println("Difference between the two inner products: ", in1 - in2)



win2 = normalize_window(win, hop)
padding=true
MTMons = z -> STFT.istft(STFT.stft(Complex.(z), win2, hop; padding=padding
            ) , win2, hop; padding=padding
            )

P = MTM(ones(signal_len)) .|> real
stft(ones(Complex, signal_len), win, hop) .|> abs2 .|> log10 |> heatmap
stft(Complex.(x), win, hop) .|> abs2 .|> log10 |> heatmap
plot(P[10:end-10])
P2 = MTMons(ones(signal_len)) .|> real
plot!(P2[10:end-10])


plot([P.*x, real.(P.*x - MTM(x))])
plot(real.(P.*x - MTM(x)))



# initialize the normalization window
norm = win.^2
for i = 1+hop:hop:win_len
    norm[1:end-i+1] += win[i:end].^2
    norm[i:end] += win[1:end-i+1].^2
end

win_hat = win ./ sqrt.(real(norm))

plot(win, label = "original window function")
plot!(win_hat, label = "modified window function")

win_hatons = STFT.normalize_window(win, hop)
win_hatons2 = STFT.normalize_window(win_hatons, hop)
plot!(win_hatons)
plot!(win_hatons2)
(win_hat - win_hatons ).|> abs2 |> sum
(win_hatons2 - win_hatons ).|> abs2 |> sum

     
# modify M2
M̂2 = Z -> Map2(Z, win_hat)
M̂2T = M̂2

# modify M
M̂ = z -> M3(M̂2(M1(z)))
M̂T = Z -> M1T(M̂2T(M3T(Z)))

N = win_len - hop
P̂ = M̂T(M̂(ones(signal_len)))
dif = M̂T(M̂(x)) - x

plot(x, label = "input signal")
plot!(real.(P̂), label = "P̂")
plot!(real(dif), label = "real part of the reconstruction error")
plot!(imag(dif), label = "imaginary part of the reconstruction error")

plot(win)
plot!(win_hatons)

win2 = win_hatons
P̂ = STFT.istft(STFT.stft(ones(Complex, signal_len), win2, hop) , win2, hop)
dif = STFT.istft(STFT.stft(Complex.(x), win2, hop) , win2, hop) - x
x2 = STFT.istft(STFT.stft(Complex.(x), win2, hop) , win2, hop)
plot(x, label = "input signal")
plot!(real.(P̂), label = "P̂")
plot!(real(dif), label = "real part of the reconstruction error")
plot!(imag(dif), label = "imaginary part of the reconstruction error")
plot!(x2 .|> real)
abs2.(x - x2) |> sum


win = hann_window(win_len)
winnorm = normalize_window(win, hop)
plot(win)
plot!(winnorm)

win2 = winnorm
padding = true

x = randn(signal_len)
plot(x)
X = rstft(x, win2, hop, padding=padding)
plot!(irstft(X, win2, hop, padding=padding))
x - irstft(X, win2, hop, padding=padding) .|> abs2 |> sum
heatmap(X .|> abs2 .|> log10)
Y = zeros(Complex, size(X))
Y = Y .+ 1e-10
Y[1:5, :] = X[1:5, :]
heatmap(Y .|> abs2 .|> log10)
y = irstft(Y, win2, hop, padding=padding)
x - y .|> abs2 |> sum
plot(x)
plot!(y)
Z = rstft(y, win2, hop, padding=padding)
X-Y .|> abs2 |> sum
X-Z .|> abs2 |> sum
Y-Z .|> abs2 |> sum
heatmap(Y .|> abs)
heatmap(Z .|> abs)
heatmap(Y-Z .|> abs)

M = abs.(X .|> abs)
heatmap(M)

m = irmstftm(M, win2, hop)
plot(x)
plot!(m)
heatmap(M)
heatmap(rstft(m, win2, hop) .|> abs)




plot((hann_window(100)))
plot!(normalize_window((hann_window(100)),39))

sum(normalize_window(hann_window(50), 15))
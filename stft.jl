module STFT

using FFTW
using Plots
using WAV

export stft, rstft, istft, irstft, imstft, irmstft, imstftm, irmstftm, normalize_window, hann_window

#inspiratie van https://github.com/s-zymon/STFT.jl

#uitleg window normalisering:
# https://github.com/ilkerbayram/STFT-Julia/blob/master/STFT_notes.ipynb

# zie ook:
# https://github.com/ilkerbayram/Basic-STFT/blob/master/GriffinLim.ipynb

# zie ook de grifin-lim paper
# Signal Estimation from Modified Short-Time Fourier Transform DANIEL W . GRIFFIN A N D J A E S. LIM, SENIOR MEMBER, IEE

# zie ook:
# https://www.audiolabs-erlangen.de/resources/MIR/FMP/C2/C2_STFT-Inverse.html
# https://www.audiolabs-erlangen.de/resources/MIR/FMP/C8/C8S1_SignalReconstruction.html


_fft(x::AbstractMatrix{<:Real}, d) = rfft(x, d)
_fft(x::AbstractMatrix{<:Complex}, d) = fft(x, d)


function maaksegments(x :: Vector{T}, window_length :: Int, hop :: Int; padding=true) where T<:Number #:: Array{<:Number, 2}
    
    signal_length = length(x)
    num_segments = div(signal_length - window_length, hop)+1

    numpadsegs = div(window_length - hop, hop, RoundUp)
    num_segments = padding ? (num_segments + 2 *numpadsegs) : num_segments
    padspul = zeros(T, hop * numpadsegs)

    #X = Array{T, 2}(undef, (window_length, num_segments))
    X = zeros(T, (window_length, num_segments))
    xalt = padding ? [padspul ; x ; padspul] : x
    for j = 1:num_segments
        startindex = (j-1)*hop+1
        eindindex = startindex + window_length - 1
        X[:, j] = xalt[startindex:eindindex]
    end
    return X
end

function applyWindow(X :: Array{<:Number, 2}, w :: Vector{<:Number}) #:: Array{Number, 2}
    return X .* w
end

"calculates unitary stft, for complex vectors calculates twosided, for real vectors calculates onesided (resulting in n/2 + 1 frequency bins)"
function stft(x :: Vector{<:Number}, w :: Vector{<:Number}, hop :: Int; padding=true) :: Array{<:Complex, 2}
    window_length = length(w)
    segs = maaksegments(x, window_length, hop; padding=padding)
    segsw = applyWindow( segs, w )
    return _fft(segsw, 1) ./ sqrt(length(w))
end

"calculates unitary stft for real vector, so onesided (resulting in n/2 + 1 frequency bins)"
function rstft(x :: Vector{<:Real}, w :: Vector{<:Number}, hop :: Int; padding=true) :: Array{<:Complex, 2}
    stft(x, w, hop; padding=padding)
end


function overlapadd(X :: Array{T, 2}, win, hop :: Int) :: Vector{T} where T <: Number
    X = X .* win

    num_segments = size(X)[2]
    window_length = size(X)[1]
    signal_length = hop * (num_segments - 1) + window_length
    x = zeros(T, signal_length) :: Vector{T}
    #d = zeros(T, signal_length) :: Vector{T}
    for j = 1:num_segments
        startindex = (j-1)*hop+1
        eindindex = startindex + window_length - 1
        x[startindex:eindindex] += X[:, j]# .* win
        #d[startindex:eindindex] += win # .^2
    end
    return x# ./ d
end


function istftcommon(X :: Array{<:Number, 2}, w :: Vector{<:Number}, hop :: Int, fftfun; padding=true) # :: Vector{<:Number}
    win_len = length(w)

    numpadsegs = div(win_len - hop, hop, RoundUp)
    #if padding 
    #    #numpadsegs = div(win_len - hop, hop, RoundUp)
    #    #paddingsegs = zeros((size(X)[1], numpadsegs))
    #    #Xalt = [paddingsegs X paddingsegs]
    #    segs = fftfun(Xalt, win_len) .* sqrt(length(w))
    #else
        segs = fftfun(X, win_len) .* sqrt(length(w))
    #end

    #segsw = applyWindow(X, window)
    x = overlapadd(segs, w, hop)
    return padding ? x[numpadsegs*hop+1:end-numpadsegs*hop] : x
end

"calculates unitary inverse stft resulting in complex (time-domain) function.
 NB: requires a normalized window for perfect reconstruction"
function istft(X :: Array{<:Number, 2}, w :: Vector{<:Number}, hop :: Int; padding=true) :: Vector{<:Complex}
    fftfun = (Y, wlen) -> ifft(Y, 1)
    return istftcommon(X, w, hop, fftfun, padding=padding)
end

"calculates unitary inverse stft resulting in real (time-domain) function.
 NB: requires a normalized window for perfect reconstruction"
function irstft(X :: Array{<:Number, 2}, w :: Vector{<:Number}, hop :: Int; padding=true) :: Vector{<:Real}
    fftfun = (Y, wlen) -> irfft(Y, wlen, 1)
    return istftcommon(X, w, hop, fftfun, padding=padding)
end




"""
calculates least squares error (complex) signal estimate from modified (unitary) stft.
effectively calculates (unitary) inverse stft resulting in complex (time-domain) function, 
but the stft may have been modified and there may not exist a signal that results in the
modified stft
"""
function imstft()
    error("unimplemented: just use regular istft with a normalized window, is same")
end


"""
calculates least squares error (real) signal estimate from modified (unitary) stft.
effectively calculates (unitary) inverse stft resulting in complex (time-domain) function, 
but the stft may have been modified and there may not exist a signal that results in the
modified stft
"""
function irmstft()
    error("unimplemented: just use regular irstft with a normalized window, is same")
end



function imstftmcommon(X :: Array{<:Real, 2}, w :: Vector{<:Number}, hop :: Int, stftfun, istftfun; padding=true) # :: Vector{<:Number}
    #heatmap(X .|> log10) |> display
    #Y = Complex.(X) #dit is gewoon met 0 als complexe component, heeft dus zelfde magnitude
    Y = cispi.(rand(Float64, size(X))) .* X

    #heatmap(Y .|> abs .|> log10) |> display
    # TODO global minimum is niet gegarandeerd. Initial guess dus beter over nadenkn
    # of misschien meerdere opties als parameter aanbieden

    
    println("start MSE $( (X - abs.(Y)) .|> abs2 |> sum)")

    for i = 1:100
        y = istftfun(Y, w, hop; padding)
        Y = stftfun(y, w, hop; padding)
        #heatmap(Y .|> abs .|> log10) |> display
        #wavplay(y, 48000)
        println("current MSE $( (X - abs.(Y)) .|> abs2 |> sum)")
        angles = angle.(Y)
        Y = angles .* X
        #heatmap(Y .|> abs .|> log10) |> display
    end
    
    y = istftfun(Y, w, hop; padding)
    return y
end

"""
calculates least squares error (complex) signal estimate from modified (unitary) stft magnitudes.
effectively calculates (unitary) inverse stft resulting in complex (time-domain) function, 
but we only have magnitudes from the stft and these may have been modified and there may not 
exist a signal that results in the modified 

NB: magnitude is `abs`! not `abs2` or logarithmt
"""
function imstftm(X :: Array{<:Real, 2}, w :: Vector{<:Number}, hop :: Int; padding=true) :: Vector{<:Complex}
    istftfun = istft
    stftfun = stft
    imstftmcommon(X, w, hop, stftfun, istftfun, padding=padding)
end


"""
calculates least squares error (real) signal estimate from modified (unitary) stft magnitudes.
effectively calculates (unitary) inverse stft resulting in real (time-domain) function, 
but we only have magnitudes from the stft and these may have been modified and there may not 
exist a signal that results in the modified 

NB: magnitude is `abs`! not `abs2` or logarithmt
"""
function irmstftm(X :: Array{<:Real, 2}, w :: Vector{<:Number}, hop :: Int; padding=true) :: Vector{<:Real}
    istftfun = irstft
    stftfun = rstft
    imstftmcommon(X, w, hop, stftfun, istftfun, padding=padding)
end




"normalize window such that the stft is perfectly representable,
scales window such that \\sum_k w^2(k*S-n) = 1
"
function normalize_window(w, hop)
    win_len = length(w)
    w2 = abs2.(w)
    P = copy(w2)
    for i = 1 + hop : hop : win_len
        P[i:end] += w2[1:end-i+1] # voor windows die naar rechts geschoven zijn
        P[1:end-i+1] += w2[i:end] # voor windows die naar links geschoven zijn
    end
    return w ./ sqrt.(P)
end


"hann window"
function hann_window(window_length)
    sin.(π * (1:window_length)/(window_length+1)).^2;
end



end # module
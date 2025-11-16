module STFT

using FFTW

export stft, rstft, istft, irstft, imstft, irmstft, imstftm, irmstftm, normalize_window, hann_window

#inspiratie van https://github.com/s-zymon/STFT.jl

#uitleg window normalisering:
# https://github.com/ilkerbayram/STFT-Julia/blob/master/STFT_notes.ipynb

# zie ook:
# https://github.com/ilkerbayram/Basic-STFT/blob/master/GriffinLim.ipynb

# zie ook de grifin-lim paper
# Signal Estimation from Modified Short-Time Fourier Transform DANIEL W . GRIFFIN A N D J A E S. LIM, SENIOR MEMBER, IEE


_fft(x::AbstractMatrix{<:Real}, d) = rfft(x, d)
_fft(x::AbstractMatrix{<:Complex}, d) = fft(x, d)


function maaksegments(x :: Vector{T}, window_length :: Int, hop :: Int) where T<:Number #:: Array{<:Number, 2}
    signal_length = length(x)
    num_segments = div(signal_length - window_length, hop)+1
    X = Array{T, 2}(undef, (window_length, num_segments))
    for j = 1:num_segments
        startindex = (j-1)*hop+1
        eindindex = startindex + window_length - 1
        X[:, j] = x[startindex:eindindex]
    end
    return X
end

function applyWindow(X :: Array{<:Number, 2}, w :: Vector{<:Number}) #:: Array{Number, 2}
    return X .* w
end

"calculates unitary stft, for complex vectors calculates twosided, for real vectors calculates onesided (resulting in n/2 + 1 frequency bins)"
function stft(x :: Vector{<:Number}, w :: Vector{<:Number}, hop :: Int) :: Array{<:Complex, 2}
    window_length = length(w)
    segs = maaksegments(x, window_length, hop)
    segsw = applyWindow( segs, w )
    return _fft(segsw, 1) ./ sqrt(length(w))
end

"calculates unitary stft for real vector, so onesided (resulting in n/2 + 1 frequency bins)"
function rstft(x :: Vector{Real}, w :: Vector{<:Number}, hop :: Int) :: Array{<:Complex, 2}
    stft(x, w, hop)
end


function overlapadd(X :: Array{T, 2}, hop :: Int) :: Vector{T} where T <: Number
    num_segments = size(X)[2]
    window_length = size(X)[1]
    signal_length = hop * (num_segments - 1) + window_length
    x = zeros(T, signal_length) :: Vector{T}
    for j = 1:num_segments
        startindex = (j-1)*hop+1
        eindindex = startindex + window_length - 1
        x[startindex:eindindex] += X[:, j]
    end
    return x
end


function istftcommon(X :: Array{<:Number, 2}, w :: Vector{<:Number}, hop :: Int) # :: Vector{<:Number}
    segsw = applyWindow(X, w)
    return overlapadd(segsw, hop)
end

"calculates unitary inverse stft resulting in complex (time-domain) function"
function istft(X :: Array{<:Number, 2}, w :: Vector{<:Number}, hop :: Int) :: Vector{<:Complex}
    segs = ifft(X, 1) .* sqrt(length(w))
    return istftcommon(segs, w, hop)
end

"calculates unitary inverse stft resulting in real (time-domain) function"
function irstft()
    error("unimplemented")
end




"""
calculates least squares error (complex) signal estimate from modified (unitary) stft.
effectively calculates (unitary) inverse stft resulting in complex (time-domain) function, 
but the stft may have been modified and there may not exist a signal that results in the
modified stft
"""
function imstft()
    error("unimplemented")
end


"""
calculates least squares error (real) signal estimate from modified (unitary) stft.
effectively calculates (unitary) inverse stft resulting in complex (time-domain) function, 
but the stft may have been modified and there may not exist a signal that results in the
modified stft
"""
function irmstft()
    error("unimplemented")
end





"""
calculates least squares error (complex) signal estimate from modified (unitary) stft magnitudes.
effectively calculates (unitary) inverse stft resulting in complex (time-domain) function, 
but we only have magnitudes from the stft and these may have been modified and there may not 
exist a signal that results in the modified 
"""
function imstftm()
    error("unimplemented")
end


"""
calculates least squares error (real) signal estimate from modified (unitary) stft magnitudes.
effectively calculates (unitary) inverse stft resulting in real (time-domain) function, 
but we only have magnitudes from the stft and these may have been modified and there may not 
exist a signal that results in the modified 
"""
function irmstftm()
    error("unimplemented")
end




"normalize window such that the stft is perfectly representable"
function normalize_window()
    error("unimplemented")
end


"hann window"
function hann_window(window_length)
    sin.(π * (1:window_length)/(window_length+1)).^2;
end



end # module
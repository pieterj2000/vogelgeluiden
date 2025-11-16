module STFT

using FFTW

export stft, rstft, istft, irstft, imstft, irmstft, imstftm, irmstftm

#inspiratie van https://github.com/s-zymon/STFT.jl

#uitleg window normalisering:
# https://github.com/ilkerbayram/STFT-Julia/blob/master/STFT_notes.ipynb

# zie ook:
# https://github.com/ilkerbayram/Basic-STFT/blob/master/GriffinLim.ipynb

# zie ook de grifin-lim paper
# Signal Estimation from Modified Short-Time Fourier Transform DANIEL W . GRIFFIN A N D J A E S. LIM, SENIOR MEMBER, IEE

"calculates unitary stft, for complex vectors calculates twosided, for real vectors calculates onesided (resulting in n/2 + 1 frequency bins)"
function stft()
    error("unimplemented")
end

"calculates unitary stft for real vector, so onesided (resulting in n/2 + 1 frequency bins)"
function rstft()
    error("unimplemented")
end

"calculates unitary inverse stft resulting in complex (time-domain) function"
function istft()
    error("unimplemented")
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






end # module
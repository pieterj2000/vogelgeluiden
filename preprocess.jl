using JSON
using StatsBase
using FFMPEG_jll
using WAV

using Plots
using STFT

files = readdir("data/downloads/")
audiofiles = filter(x -> in('.', x), files)
metafiles = filter(x -> !in('.', x), files)
uitgangen = unique(map(f -> split(f, '.')[end], audiofiles))

samplerates = sort(map( x-> JSON.parsefile("data/downloads/" * x)["smp"],  metafiles)) 

countmap(samplerates)
# er zijn 472 dingen met 48k rate, dus die pakken we. 
# TODO: we zouden ook nog heel makkelijk de 96k kunnen pakken en resapmlen (of zelfs in twee samples veranderen door de even samples en de oneven samples te pakne)
# TODO kijken of lagere sample rate doebaar is. Dan zouden we met 24k er 5 samples bij krijgen, of met 16k nog 2 samples

samplerate = 48000 :: Int

goedeaudiofiles = filter( x -> JSON.parsefile("data/downloads/" * split(x, ".")[1])["smp"] == string(samplerate), audiofiles)

file = goedeaudiofiles[1]

"sdfsdlfjsdlfkj"
function readwav(file :: String) 
    filepath, fileio = mktemp()
    run(`$(ffmpeg()) -i $file -f wav $(filepath) -y`)
    #println(filepath)
    wav, rate = wavread(filepath)
    rm(filepath)
    return wav
end


# x = load("data/downloads/XC1049551.wav")
# x = load("data/downloads/XC1004643.mpeg")
# x = load("data/downloads/XC1052236.mp3")

x = readwav("data/downloads/XC1004643.mpeg")


y = x[1:(10*samplerate),1]

wavplay(y, samplerate)




win = win_len -> sin.(π * (1:win_len)/(win_len+1)).^2
winnorm = win_len -> begin
        win1 = win(win_len)
        norm = win1.^2
        for i = 1+hop:hop:win_len
            norm[1:end-i+1] += win1[i:end].^2
            norm[i:end] += win1[1:end-i+1].^2
        end

        win_hat = win1 ./ sqrt.(real(norm))  
        return win_hat
    end




win_len = 50
win = sin.(π * (1:win_len)/(win_len+1)).^2;
hop = 15
plot(1:win_len, win)



W = 4800
#w = ones(W)
H = 600
L = W - H

Y = stft(y, winnorm(W), L)
s = log10.(abs2.(Y)) 
s = s .- maximum(s)
s = clamp.(s, -5,0)
heatmap(s, color=cgrad(:grays, rev=true))

Y2 = 10 .^ Y

y2 = istft(Y, winnorm(W), L)

plot([y[1:1000] y2[1:1000]])

wavplay(y, samplerate)
wavplay(y2, samplerate)


4test = ones(50+29*15)
Test = stft(test, winnorm(50), 50-15)


heatmap(abs2.(Test) .|> log10)
testr = istft(Test, winnorm(50), 50-15)
plot((testr))

maximum(abs.(testr))

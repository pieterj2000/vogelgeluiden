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

x[1:1000, :]

wavplay(x[1:100000, :], samplerate)


y = x[:,1]


W = 64
w = ones(W)
H = 10
L = W - H

Y = stft(y, w, L)
s = abs2.(Y)
heatmap(s)
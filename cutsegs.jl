using JSON
using StatsBase
using FFMPEG_jll
using WAV
using Plots



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

function readwav(file :: String) 
    filepath, fileio = mktemp()
    run(`$(ffmpeg()) -i $file -f wav $(filepath) -y`)
    #println(filepath)
    wav, rate = wavread(filepath)
    close(fileio)
    rm(filepath; force=true)
    return wav
end

function conv(f, g)
    glen = length(g)
    glen2 = div(glen,2)
    pad = zeros(glen2)
    fp = [pad ; f ; pad]
    result = zeros(length(f))
    for i=1:length(f)
        result[i] = fp[i:(i+glen-1)]' * g
    end
    return result
end
function conv(f, g)
    glen = length(g)
    glen2 = div(glen,2)
    pad = zeros(glen2)
    fp = [pad ; f ; pad]
    result = zeros(length(f))
    for i=1:length(f)
        result[i] = fp[i:(i+glen-1)]' * g
    end
    return result
end


file = rand(goedeaudiofiles)#[1]
x = readwav("data/downloads/$file")
plot(x)
y = x[:, 2]
mask = abs.(y) .> 0.03
plot(y)
plot!(mask)


plot(y[1:48000])
plot!(mask[1:48000])

wavplay(y[1:48000], samplerate)


plot(mask[11800:11900])
mask[11830:29000]
maskstarts[11830:29000]

wavplay(y, samplerate)
maskstarts = conv(mask, [-1,1]) .|> Int



mingat = 50
indices = findall(i -> i != 0, maskstarts)
vorigestart = indices[1]
vorigeeind = indices[1]
segs = []
#its = Iterators.takewhile(<(48000), Iterators.drop(indices, 1))
its = Iterators.drop(indices, 1)
for i in its
    if maskstarts[i] == 1
        # start
        if (i - vorigeeind) > 48*mingat  ## groot gat
            push!(segs, (vorigestart-mingat*24, vorigeeind+mingat*24))
            vorigestart = i
        end # anders hoeven we niets te doen...
    else 
        # eind
        vorigeeind = i
    end
end
push!(segs, (vorigestart-mingat*24, vorigeeind+mingat*24))
segs

plot(mask)
plot!(y)
#plot!(mask[1:48000])
map( t-> plot!(t[1]:t[2], x -> rand() .* 0.4), segs)
plot!([0]; legend=false)


convmask = conv(max.(y, mask), ones(48*10)) ./ (48*10)
plot!(convmask)



numfiles = length(goedeaudiofiles)
sigs = Vector{Float64}[]
fours = Array{ComplexF64,2}[]
for i=1:numfiles
    println(i)
    file = goedeaudiofiles[i]
    x = readwav("data/downloads/$file")
    if length(x[:,1]) < 480000
        continue
    end
    y = x[1:(10*samplerate),1]
    Y = rstft(y, win, H)
    push!(sigs, y)
    push!(fours, Y)
end
mean(sigs)
plot(mean(sigs))
heatmap(mean(fours))

mean(sigs) - meansig

meansig = zeros(480000)
meanfour = zeros((601,1603))

meansig ./ 472
meanfour ./ 472


plot(meansig)
wavplay(meansig ./ 5, samplerate)
heatmap(meanfour .|> abs2 .|> log10)
heatmap(rstft(meansig, win, H) .|> abs2 .|> log10)
m = mean(meanfour .|> abs2 .|> log10)
meanfour .|> abs2 .|> log10
heatmap((meanfour .|> abs2 .|> log10) .- m)



# x = load("data/downloads/XC1049551.wav")
# x = load("data/downloads/XC1004643.mpeg")
# x = load("data/downloads/XC1052236.mp3")

#x = readwav("data/downloads/XC1004643.mpeg")
#x = readwav("piano2.wav")
file = rand(goedeaudiofiles)
#file = goedeaudiofiles[6]
x = readwav("data/downloads/$file")


y = x[1:(10*samplerate),1]
#y = x #x[1:(3*samplerate),1]



W = 1200
#w = ones(W)
H = 300
L = W - H


win = hann_window(W)
win = normalize_window(win, H)
# plot(hann_window(W))
# plot!(win)

Y = rstft(y, win, H)
#Y = rstft(x[:,1], win, H)

s = abs.(Y)
s = log10.(s)
#s = s .- maximum(s)
#s = clamp.(s, -5,0)
heatmap(s, color=cgrad(:grays, rev=true))
#heatmap(Y .|> abs2 .|> log10, color=cgrad(:grays, rev=true))

powers = ones(size(Y)[1])' * abs2.(Y)
scale = 600 / mean(powers) / 20
scale2 = 600 / maximum(powers)
#plot!(powers'*scale) |> display
plot!(powers'*scale2) |> display
powersconv = conv(powers'*scale2, ones(9)./9)
plot!(powersconv)
#plot!((1:480000) ./ 480000 .* 1603, y .*400 .- 200)
yloud = (abs.(y) .> 0.03) .* y
plot!((1:480000) ./ 480000 .* 1603, y .*400 .- 200)
plot!((1:480000) ./ 480000 .* 1603, yloud .*400 .- 200)
wavplay(y, samplerate)
file

Yvals = reshape(abs2.(Y) .|> log10, prod(size(Y)))
ec = ecdf(Yvals)
r = sort(Yvals)[begin:10:end]
rvals = map(ec, r)
#r = range(minimum(Yvals), maximum(Yvals),100)
plot(r,rvals)
quantile(r, 0.75)


plot(yloud)
y2 = conv(yloud, ones(100*div(samplerate, 1000))) ./ 50
plot!(y2)



Yvals = (powers' * scale2)
ec = ecdf(Yvals)
r = sort(Yvals)[begin:end]
rvals = map(ec, r)
#r = range(minimum(Yvals), maximum(Yvals),100)
plot(r,rvals, xscale=:log10)
quantile(r, 0.75)


function conv(f, g)
    glen = length(g)
    glen2 = div(glen,2)
    pad = zeros(glen2)
    fp = [pad ; f ; pad]
    result = zeros(length(f))
    for i=1:length(f)
        result[i] = fp[i:(i+glen-1)]' * g
    end
    return result
end
#conv(ones(10), ones(5)./5)

powersconv = conv(powers'*scale, ones(9)./9)
plot!(powersconv)


wavplay(y, samplerate)

wavplay(y[1:(3*samplerate), 1], samplerate)



plot(y .+ 100)


wavplay(y, samplerate)


YCrop = ones(Complex, size(Y)) * 1e-20
YCrop[150:end, :] = Y[150:end, :]
ycrop = irstft(YCrop, win, H)
plot(y)
plot!(ycrop)
heatmap(abs2.(Y) .|> log10)
heatmap(abs2.(YCrop) .|> log10)
heatmap(abs2.(rstft(ycrop, win, H)) .|> log10)
wavplay(y, samplerate)
wavplay(ycrop, samplerate)



Y2 = 10 .^ Y

y2 = irstft(Y, win, H)

abs2.(x[1:1720200, 1] - y2) |> sum
abs2.(x[4800:1720200-4800, 1] - y2[4800:1720200-4800]) |> sum

plot([x[1:10000, 1] y2[1:10000]])
x[1:10000, 1] - y2[1:10000] .|> abs2 |> sum


M = abs.(Y)
heatmap(M .|> log10, color=cgrad(:grays, rev=true))
m = irmstftm(M, win, H)

M2 = Complex.(M)
M2 = cis.(angle.(Y)) .* M
M2 = cispi.(rand(Float64, size(M))) .* M
mslecht = irstft(M2, win, H)
plot(y)
plot!(mslecht)
plot!(m)
wavplay(y, samplerate)
wavplay(mslecht, samplerate)
wavplay(m, samplerate)





sum(abs2.(y))
sum( abs2.(stft(Complex.(y), win, H)))


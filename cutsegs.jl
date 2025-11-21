using JSON
using StatsBase
using FFMPEG_jll
using WAV
using Plots
using DSP



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
    run(`$(ffmpeg()) -i $file -f wav $(filepath) -y`, devnull, devnull, devnull)
    #println(filepath)
    wav, rate = wavread(filepath)
    close(fileio)
    rm(filepath; force=true)
    return wav
end


## Deze vervangen we met DSP.jl conv
# function conv(f, g)
#     glen = length(g) :: Int
#     glen2 = div(glen,2)
#     pad = zeros(glen2)
#     fp = [pad ; f ; pad]
#     result = zeros(length(f))
#     for i=1:length(f)
#         result[i] = fp[i:(i+glen-1)]' * g
#     end
#     return result
# end

function getsegs(y :: AbstractVector) 
    #y = y[1:200000]

    #mask = abs.(y) .> 0.06  #0.03
    #plot(y)
    #plot!(mask)
    yalt = conv(abs.(y), ones(751)./ 50)
    pad = div(751,2)
    yalt = yalt[1+pad:end-pad] # dit hoeft alleen maar omdat de conv van DSP zero-pad
    thresh = 0.10 + quantile(yalt, 0.33)*4
    maskalt = yalt .> thresh #0.2  #0.03
    #quantile(yalt, 0.33)
    #median(yalt)

    # #y[170000:end-30000]
    # #plot(y[179880:end-60050])
    # #plot!(mask[179880:end-60050])
    # #plot!(yalt[179880:end-60050])
    # println(size(y))
    # plot(y) #[130000:135000])
    # #plot!(mask)
    # plot!(maskalt .* rand(length(y)) .* 0.4 .+ 0.2)
    # plot!(yalt) |> display


    #maskalt
    maskstarts = conv(maskalt, reverse([-1,1,0])) .|> Int 
    # wij deden eigenlijk cross correlation, daarom is die reverse hierboven nodig
    pad = div(3,2)
    maskstarts = maskstarts[1+pad:end-pad] # dit hoeft alleen maar omdat de conv van DSP zero-pad
    # 50 was al enigszins prima goeie,
    # TODO: misschien een tweede laag over doen met langere tijd
    mingat = 250  # 250 werkte ook goed, maar bij heel goede opnames, te goed
    indices = findall(i -> i != 0, maskstarts)
    if isempty(indices)
        println("Geen segmenten gespot")
        return []
    end

    vorigestart = indices[1]
    vorigeeind = indices[1]
    segs = []
    #its = Iterators.takewhile(<(48000), Iterators.drop(indices, 1))
    its = Iterators.drop(indices, 1)
    for i in its
        #println("$i: $(maskstarts[i])")
        if maskstarts[i] == 1
            # start
            if (i - vorigeeind) > 48*mingat  ## groot gat
                push!(segs, (vorigestart-mingat*24, vorigeeind+mingat*24))
                #push!(segs, (vorigestart, vorigeeind))
                vorigestart = i
            end # anders hoeven we niets te doen...
        else 
            # eind
            vorigeeind = i
        end
    end
    if maskalt[end] == 1
        vorigeeind = length(maskalt)
    end
    push!(segs, (vorigestart-mingat*24, vorigeeind+mingat*24))
    #push!(segs, (vorigestart, vorigeeind))

    # #println("klaar om te printen")
    # println(thresh)
    # println(quantile(yalt, 0.33))
    # plot([0]; legend=false)
    # println("mask geprint")
    # #plot!(mask[1:48000])
    # plot!(mask .* rand(length(y)) .* 0.4 .+ 0.2)
    # map( t-> plot!(t[1]:t[2], x -> rand() .* 0.4), segs)
    # plot!(min.(yalt, 0.5), color=:black,legend=false)
    # plot!(y, color=:blue,legend=false) |> display
    
    return segs
end

function getsegs(y :: AbstractArray{T,2}) where T
    seglist = map(getsegs, eachcol(y))
    return vcat(seglist...)
end

function seglength((l,r))
    return r - l + 1
end


# probleembestanden
#"XC121684.mp3" # veel noise
#"XC362777.mp3" # zacht
#"XC702690.mp3" # zacht
#"XC694031.mp3" # zacht, maar super wacky


rm("data/segments"; force=true, recursive=true)
mkdir("data/segments")

function writeseg(seg, stream, file, channel, segnum)
    (filename, extension) = split(file, ".")
    l = length(stream)
    start = clamp(seg[1], 1, l)
    eind = clamp(seg[2], 1, l)
    x = stream[start:eind]
    wavwrite(x, samplerate, "data/segments/$filename-$channel-$segnum.wav")
end


allesegs = []
for (i,f) in enumerate(goedeaudiofiles)
    println("$i: $f")
    x = readwav("data/downloads/$f")
    for (channel, y) in enumerate(eachcol(x))
        segs = getsegs(x)
        wr = (seg,i) -> writeseg(seg, y, f, channel, i)
        map(wr, segs, 1:length(segs))
        println(length(segs))
        append!(allesegs, segs)
    end
end
allesegs


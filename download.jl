
import HTTP
import JSON
using RateLimiter

const key = read("key", String)

function getPage(query :: String, page :: Integer )
    spul = "https://xeno-canto.org/api/3/recordings?query="
    spulen = spul * HTTP.escapeuri(query) * "&key=" * key * "&per_page=500&page=" * string(page)
    #println(spulen)
    r = HTTP.get(spulen)
    b = JSON.parse(String(r.body))

    numpages = b["numPages"]
    recs = b["recordings"]
        
    if page == numpages
        return recs
    else
        return [ recs ; getPage(query, page+1) ]
    end
end

function getQuery(query :: String)
    getPage(query, 1)
end

function makedownloadsfolder()
    rm("data/downloads"; force=true, recursive=true)
    mkpath("data/downloads")
end


function metaToFile(ding :: JSON.Object{String, Any})
    id = ding["id"] # is opgeslagen als string
    open("data/downloads/XC" * id, "w") do file
        JSON.json(file, ding, pretty=true)
    end
end

function metaToFile(spul :: Vector)
    metaToFile.(spul)
end

limiter = TokenBucketRateLimiter(10,3,0)

function download(ding :: JSON.Object{String, Any})
    url = ding["file"]
    r = @rate_limit limiter 1 HTTP.get(url, status_exception=false)
    if r.status != 200 
        println("ERROROROROROR status code " * string(r.status) * " bij url: " * url)
        rm("data/downloads/XC" * ding["id"])
    else
        filenamecont = HTTP.header(r, "content-disposition")
        filenameparts = split(filenamecont, ['.', ' ', '\"'])
        filename = filenameparts[3] * "." * filenameparts[end-1]
        #if filenameparts[end-1] ∉ ["wav", "mp3"]
        #    println("################################## " * filename)
        #end
        println(filename)
        write("data/downloads/" * filename, HTTP.body(r))
    end
end

dingen = getQuery("sp:\"parus major\" q:A")

makedownloadsfolder()
metaToFile(dingen)
download.(dingen)


exit()
## als je los wilt proberen
dingen = getQuery("nr:125484")
metaToFile.(dingen)
download.(dingen)

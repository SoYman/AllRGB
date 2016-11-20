println("Script started!")

using Images
using Images.data
using ImageView
using Color
using FixedPointNumbers
# using GLPlot, GLAbstraction, GLFW
using Debug
using LowDimNearestNeighbors
using LowDimNearestNeighbors.Result
using Distances
using ProfileView
using IProfile


# Constants
bitdepth = 6
winy = 1 << int(bitdepth*1.5-0.5)
winx = 1 << int(bitdepth*1.5)
bitstride = uint8(256 / (1 << bitdepth))

# Utility functions
luminance(color) = convert(LAB, color).l
hue(color) = convert(HSV, color).h
value = hue

global running = true

<(a::RGB, b::RGB) = value(a) < value(b)
>(a::RGB, b::RGB) = value(a) > value(b)
# ==((a1::Uint8, a2::Uint8, a3::Uint8), (b1::Uint8, b2::Uint8, b3::Uint8)) = a1 == b1 && a2 == b2 && a3 == b3


# Unordered code
println("Constants initialised")
# colors = [e for e in [(uint8(r)::Uint8,uint8(g)::Uint8,uint8(b)::Uint8) for r = 0:bitstride:255, g = 0:bitstride:255, b = 0:bitstride:255]]
# println("Colors array initialised")
# preprocess!(colors)
# println("Colors array pre processed")

# sort!(colors, alg=QuickSort, lt=<) # Do you even sort brah!?
# println("Colors array sorted")
# window = createdisplay(;async=true, w=winx, h=winy)
# println("Window created")
colored = [false for x in 1:winx, y in 1:winy]
locdone = falses(winx, winy)
# coldone = falses(bitdepth, bitdepth, bitdepth)
# result = Texture(Vec4[Vec4(0,0,0,1)for i=1:winx, j=1:winy])
finalimg = Image([RGB{Ufixed8}(0,0,0) for x = 1:winx, y = 1:winy])
todraw = [(div(winx,2), div(winy,2))]
# col = (uint8(255/bitstride),uint8(255/bitstride),uint8(255/bitstride))
col = (0,0,0)

offsets1 = [(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,1),(1,-1),(-1,-1)]
offsets2 = [(1,-2),(-1,2),(2,1),(-2,-1),(-1,-2),(1,2),(2,-1),(-2,1),(-1,-1),(1,1),(1,-1),(-1,1)]

function getindex(b::BitArray{2}, c::(Int,Int))
    return b[c[1],c[2]]
end

@debug function grow(locations)
    result = Array((Int64,Int64),0)
    # @bp
    if length(locations) == 0
        locations = [ind2sub(size(finalimg),findfirst(locdone))]
    end
    for l in locations
        for o in offsets1
            c = (l[1] + o[1], l[2] + o[2])
            if locdone[c]
                continue
            else
                push!(result, c)
            end
        end
    end
    return result
end

function getcolor(current)
    t = hindex(bitdepth, current)
    # print(t," ")
    res = hindexInverse(bitdepth, t+1)
    if sum(res)-sum(current) > 1
        println(res, " ", current, " ", t)
    end
end

function iterate(todraw, color)
    c = color
    temp = randsubseq(todraw, 0.6)
    if length(temp) == 0
        temp = todraw
    end
    for (x,y) in temp
        locdone[x,y] = true
        c = getcolor(color)
        # pop!(colors,(r,g,b))
        # deleteat!(colors, c)
        finalimg[x,y] = RGB{Ufixed8}(c[1]*bitstride,c[2]*bitstride,c[3]*bitstride)
    end
    todraw = setdiff(todraw, temp)
    todraw = grow(todraw)
    return  finalimg, todraw, c
end



closeView(event) = global running = false

using Debug
using ProfileView
# Hilbert logic from https://www.cs.dal.ca/sites/default/files/CS-2006-07.pdf
# Algorithm a direct implementation of the logic presented there.
# returns the index into the rgb cube indexed by the 3d hilbert curve
function hindex(m, p)
    dim = 3
    h = 0
    e = 0
    d = 0
    w,l = 0,0
    for i in reverse(0:m-1)
        l = bit(p[3], i) << 2 | bit(p[2], i) << 1 | bit(p[1], i)
        l = cycleLeft(l, 1, dim)
        l = T(e, d, dim, l)
        w = grayCodeInverse(dim, l)
        e = e $ (cycleLeft(E(w), d+1, dim))
        d = (d + D(w, dim) + 1) % dim
        h = (h << dim) | w
  end
    return h
end

function hindexInverse(m, h)
    dim = 3
    e = 0
    d = 0
    p1, p2, p3 = 0,0,0
    w,l,e,d = 0,0,0,0
    for i in reverse(0:m-1)
        w = bit(h, 3*i+2) << 2 | bit(h, 3*i+1) << 1 | bit(h, 3*i)
        l = grayCode(w)
        l = TInverse(e, d, dim, l)
        p1 |= bit(l, 0) << i
        p2 |= bit(l, 1) << i
        p3 |= bit(l, 2) << i
        e = e $ (cycleLeft(E(w), d+1, dim))
        d = (d + D(w, dim) + 1) % dim
    end
    return (p1,p2,p3)
end

# trailing set bits
function tsb(i)
    r = 0
    # Split on last 4 bits
    while ((i & 7) == 7)
        # All 4 set
        r += 4 
        i >>= 4
    end
    while i & 1 == 1
        r += 1
        i >>=1
    end
    # tsb4 = 
    return r
end

function grayCodeInverse(m, g)
    i = g
    j = 1
    while j < m
        i $= (g >> j)
        j += 1
    end
    return i
end

D(i, n) = i == 0 ? 0 : ((i & 1 != 0) ? (tsb(i) % n) : (tsb(i-1) % n))

E(i) = i == 0 ? 0 : grayCode((i-1) & (~1))

grayCode(i) = i $ (i >> 1)

cycle(a, b, n) = (a >> b) | ((a << (n-b)) & (2^n - 1))

cycleLeft(a, b, n) = cycle(a, n-b, n)

T(e, d, n, b) = cycle(b$e,d+1,n)

TInverse(e, d, n, b) = T(cycle(e, d+1, n), n-d-1, n, b)

bit(n, i) = (n >> i) & 1

bitSet(n, i, val) = n | val << i

# function makeimage()
    finalimg["pixelspacing"] = [1,1]
    finalimg["spatialorder"] = ["y","x"]

    imgc, _ = view(finalimg)
    win = toplevel(imgc)
    bind(win, "<Key>", closeView)
    bind(win, "<Destroy>", closeView)

    # Profile.clear()
    # Profile.init()

    t = col
    i=1
    times = Array(Float64,1)
    while running
        tic()

        # finalimg, todraw, col = iterate(todraw, col)
        for x in 1:winx
            t = getcolor(col)
            finalimg[i,x] = RGB{Ufixed8}(col[1]*bitstride,col[2]*bitstride,col[3]*bitstride)
            col = t
        end

        view(imgc, finalimg)
        push!(times,toq())
        # sleep(0.5)
        # print(i," ")
        gc()
        if i%100 == 0
            imwrite(finalimg, "crawlimg$(int(time())).png")
        end
        i+=1
    end
    shift!(times)
    print("\nMeðal tíminn var: ", sum(times)/length(times))
    destroy(win)
    # Profile.print()
    # ProfileView.view()
    imwrite(finalimg, "crawlimg$(int(time())).png")
    # gc_enable()
    # gc()
# end
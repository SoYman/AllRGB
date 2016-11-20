println("Script started!")
# threads = 2
# nprocs() < 2 ? addprocs(2 - nprocs()) :
# print("You've got $threads processes!\n")
using Images
print("Images")
using ImageView
print(", ImageView")
using Colors
print(", Colors")
using FixedPointNumbers
print(", FixedPointNumbers")
using Debug
print(", Debug")
using Distances
print(", Distances")
using ProfileView
print(", ProfileView")
include("Blur.jl")
print(", Blur")
include("Utility.jl")
print(", Utility.\n")
print("libraries loaded\n")

# Constants

offsetsshortrook = [(1,0),(0,1)]
# offsetsking = [(1,0),(1,-1),(0,-1),(-1,-1),(-1,0),(-1,1),(0,1),(1,1)]
offsetsking = [(1,0),(1,1),(0,1),(-1,0)]

offsetsknight = [(1,2), (2,1), (2,-1), (1,-2)]

offsetsrook = ([[(n,0),(0,n)] for n in 1:8])

offsetsrook = ([[(n,0),(n,n),(0,n),(-n,0)] for n in 1:8])

offsetsten = [(y,x) for y=-10:10, x=-10:10]

# Functions

function shuffledallrgb(bitdepth)
    winy = 1 << int(bitdepth*1.5-0.5)
    winx = 1 << int(bitdepth*1.5)
    bitstride = (1 / (1 << bitdepth))


    colors = [e for e in [RGB{Float32}(r,g,b) for r = 0:bitstride:1, g = 0:bitstride:1, b = 0:bitstride:1]]
    println("Colors array initialised")
    # sort!(colors, alg=QuickSort, lt=morehue) # Do you even sort brah!?
    shuffle!(colors)
    println("Colors array permuted")
    result = Image([pop!(colors)::RGB{Float32} for y = 1:winy, x = 1:winx])


    # result = imread("img1439187330-i1930.png")


    result["pixelspacing"] = [1,1]
    result["spatialorder"] = ["y","x"]
    return result
end

shar = shuffledallrgb

function sortedallrgb(bitdepth)
    winy = 1 << int(bitdepth*1.5-0.5)
    winx = 1 << int(bitdepth*1.5)
    bitstride = uint8(256 / (1 << bitdepth))


    colors = [e for e in [RGB{Ufixed8}(r,g,b) for r = 0:bitstride:255, g = 0:bitstride:255, b = 0:bitstride:255]]
    println("Colors array initialised")
    shuffle!(colors)
    sort!(colors, alg=QuickSort, lt=morehue) # Do you even sort brah!?
    println("Colors array permuted")
    result = Image([pop!(colors)::RGB{Ufixed8} for y = 1:winy, x = 1:winx])


    # result = imread("img1439187330-i1930.png")


    result["pixelspacing"] = [1,1]
    result["spatialorder"] = ["y","x"]
    return result
end

function iterate(image, blur, sigma)

    tic()
    # shift = offsetsking[mod1(i,end)]
    # shift = offsetsten[mod1(i,end)]


    # parr1 = imfilter(image,   imagefilter, "circular")
    # parr1 = imfilter_fft(image,   gaussian2d(sigma), "circular")
    # parr1 = imfilter_LoG(image, sigma, "circular")
    # parr1 = imfilter_gaussian(image, [sigma,sigma])

    parr1 = blur(image, sigma, 2)
    offset = (0,0)
    # bestoffset = offset
    # # offsets = Set{(Int,Int)}([cat(1,offsetsking, offsetsknight)])
    # # offsets = Set{(Int,Int)}(offsetsking)
    # offsets = Set{(Int,Int)}()
    # i = 0
    # while (length(offsets) < 1) & (i < 1000)
    #     offset = (0,0)
    #     while offset == (0,0)
    #         y = randn()*sigma
    #         x = randn()*sigma
    #         # offset = (int(y < 0 ? -(-y)^0.5 : y^0.5), int(x < 0 ? -(-x)^0.5 : x^0.5))
    #         offset = (int(y/2), int(x/2))
    #     end
    #     push!(offsets, offset)
    #     offsets = setdiff(offsets, Set([-e for e in offsets]))
    #     i += 1
    # end

    farr = zeros(Float32,size(image))
    maxvariance = 0
    i=0
    # for offset in offsets
            # offset = (0,0)

        while offset == (0,0)
            y = randn()*sigma
            x = randn()*sigma
            compression = 0.75
            offset = (int(y < 0 ? -(-y)^compression : y^compression), int(x < 0 ? -(-x)^compression : x^compression))
            offset = (offset[1]%size(image)[1], offset[2]%size(image)[2])
            # offset = (int(y/2), int(x/2))
        end

        iarr2 = circshift(image, offset)
        parr2 = circshift(parr1, offset)

        farr1 = map(compare, parr1, image)
        farr2 = map(compare, parr2, iarr2)
        farr3 = map(compare, parr1, iarr2)
        farr4 = map(compare, parr2, image)
        farr  = (farr1 + farr2) - (farr3 + farr4)
        # print("$(sum(farr))  $(sum(farr))")
        maxvariance = sum(x -> x < 0 ? 0 : x, farr)
        # if maxvariance < newvariance
            bestoffset = offset
            # maxvariance = newvariance
        # end
        # i += 1
    # end

    barr = farr .> 0
    changepercent = countnz(barr)/length(barr)*100

    # final = copy(image)
    if changepercent < 80
        swapby!(image, barr, -bestoffset)
    else
        barr2 = circshift(barr, bestoffset)
        barr = barr & ~barr2
        barr2 = circshift(barr, -bestoffset)
        iarr3 = circshift(image, -bestoffset)
        image = map(elementswap, barr, iarr2, image)
        image = map(elementswap, barr2, iarr3, image)
    end
    return image, toq(), changepercent, bestoffset, maxvariance
end

function run(image, sigma = 1000)

    global running = true
    winy, winx = size(image)

    iarr1 = copy(image)

    if iarr1["spatialorder"] != ["y","x"]
        iarr1["spatialorder"] = ["y","x"]
    end


    colortype = LAB
    if colortype == LMS
        blur = gausslms
    elseif colortype == RGB
        blur = gaussrgb1d
    elseif colortype == LAB
        blur = gausslab
    elseif colortype == LUV
        blur = gaussluv
    end

    iarr1 = convert(Image{colortype{Float32}}, iarr1)
    # iarr2 = Array(colortype{Float32},(winy,winx))
    # parr1 = Array(colortype{Float32},(winy,winx))
    # parr2 = Array(colortype{Float32},(winy,winx))
    # farr  = Array(Float32,(winy,winx))
    # barr = Array(Bool,(winy,winx))

    imgc, _ = view(iarr1)
    win = toplevel(imgc)
    bind(win, "<Key>", close)
    bind(win, "<Destroy>", close)

    increment = 1
    nochange = ones(10) * 1000
    times = Array(Float32,1)
    while running
        iarr1, timeelapsed, change, offset, variance = iterate(iarr1, blur, sigma)

        pop!(unshift!(nochange, change))

        if sum(nochange) / length(nochange) < 0.1
            println("\nThere was not enough change with sigma: $sigma. ")
            sigma *= 0.9
            if sigma < 0.8
                break
            end
            nochange = ones(10) * 1000
        end

        push!(times, timeelapsed)


        Tk.update()
        if running
        view(imgc, iarr1)
        end

        print_with_color(:blue, "$increment ")
        @printf("%0.3f ", change)
        @printf("%1.2f ", variance)
        print_with_color(:green, "$offset \t")
        # gc()
        if increment%5 == 0
            println("")
        end

        if increment%10 == 0
            imwrite(iarr1, "img$(int(time()))-i$(increment).png")
        end
        increment+=1
    end

    shift!(times)
    print("\nMeðal tíminn var: ", sum(times)/length(times))
    destroy(win)

    imwrite(iarr1, "img$(int(time()))-i$(increment).png")

    allrgb(iarr1)
    println("\n", length(Set(iarr1)))
    println(length(iarr1))
    return iarr1
end

print_with_color(:green, "Done!")
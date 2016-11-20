# Inline operators

(*)(a::RGB, b::RGB) = RGB(a.r*b.r, a.b*b.b, a.g*b.g)
(-)(a::LMS, b::LMS) = LMS(a.l-b.l, a.m-b.m, a.s-b.s)
(/)(a::LMS, b::Int) = LMS(a.l / b, a.m / b, a.s / b)
(./)(a::LMS, b::Int) = (/)(a,b)

(+)(a::LAB, b::LAB) = LAB(a.l+b.l, a.a+b.a, a.b+b.b)
(.+)(a::LAB, b::LAB) = (+)(a,b)
(-)(a::LAB, b::LAB) = LAB(a.l-b.l, a.a-b.a, a.b-b.b)
(/)(a::LAB, b::Int) = LAB(a.l / b, a.a / b, a.b / b)
(./)(a::LAB, b::Int) = (/)(a,b)

(+)(a::LUV, b::LUV) = LUV(a.l+b.l, a.u+b.u, a.v+b.v)
(.+)(a::LUV, b::LUV) = (+)(a,b)
(-)(a::LUV, b::LUV) = LUV(a.l-b.l, a.u-b.u, a.v-b.v)
(/)(a::LUV, b::Int) = LUV(a.l / b, a.u / b, a.v / b)
(./)(a::LUV, b::Int) = (/)(a,b)
-(t::(Real, Real)) = (-t[1], -t[2])

# Functions

flatten{T}(a::Array{T,1}) = any(map(x->isa(x,Array),a))? flatten(vcat(map(flatten,a)...)): a
flatten{T}(a::Array{T}) = reshape(a,prod(size(a)))
flatten(a)=a

compare1(f::RGB{Float32},a::RGB{Ufixed8}) = colordiff(f, a)

compare2(f::RGB{Float32}, a::RGB{Ufixed8}) = abs(f-a)

# function compareSqEuclidean(f::RGB, a::RGB)
#     c1 = f - a
#     return (c1.r^2+c1.b^2+c1.g^2)
# end

compareEuclidean(a::RGB, b::RGB) = euclidean([a.r,a.g,a.b],[b.r,b.g,b.b])

function compareEuclidean(f::LMS, a::LMS)
    c1 = f - a
    return (c1.l^2+c1.m^2+c1.s^2)
end


compareSqEuclidean(a::LAB, b::LAB) = sqeuclidean([a.l,a.a,a.b],[b.l,b.a,b.b])

# compareEuclidean(a::LAB, b::LAB) = euclidean([a.l,a.a,a.b],[b.l,b.a,b.b])

function compareEuclidean(f::LAB, a::LAB)
    c = f - a
    return sqrt(c.l^2+c.a^2+c.b^2)
end

function compareHamming(a::RGB,b::RGB)
    return count_ones(int(a.r) $ int(b.r)) + count_ones(int(a.g) $ int(b.g)) + count_ones(int(a.b) $ int(b.b))
end

function compareMinkowski(a::RGB, b::RGB)
    minkowski([a.r,a.g,a.b],[b.r,b.g,b.b],0.5)
end

function compare5(f::RGB{Float32},a1::RGB{Ufixed8},a2::RGB{Ufixed8})
    hue = convert(LCHab, f).h
    abs(abs(hue - convert(LCHab, a1).h)-180) - abs(abs(hue - convert(LCHab, a2).h)-180)
end


morehue(a, b) = convert(HSV, a).h < convert(HSV, b).h

moreLuv(a, b) = convert(LUV, a).l < convert(LUV, b).l
morelUv(a, b) = convert(LUV, a).u < convert(LUV, b).u
moreluV(a, b) = convert(LUV, a).v < convert(LUV, b).v
moreLch(a::LCHab, b::LCHab) = (a.l < b.l)
morelCh(a::LCHab, b::LCHab) = (a.c < b.c)
morelcH(a::LCHab, b::LCHab) = (a.h < b.h)
moreLch(a, b) = convert(LCHab, a).l < convert(LCHab, b).l
morelCh(a, b) = convert(LCHab, a).c < convert(LCHab, b).c
morelcH(a, b) = convert(LCHab, a).h < convert(LCHab, b).h
morer(a, b) = a.r < b.r
moreg(a, b) = a.g < b.g
moreb(a, b) = a.b < b.b
morerg(a, b) = a.r+a.g < b.r+b.g
moregb(a, b) = a.g+a.b < b.g+b.b
morebr(a, b) = a.b+a.r < b.b+b.r
moreabs(a, b) = a.r+a.g+a.b < b.r+b.g+b.b

compare = compareEuclidean

elementswap(bool, a, b) = bool ? a : b

close(event) = global running = false

function allrgb(image)
    if length(Set(image)) == length(image) #Assert that the image only has unique colors
        print("\n\nValid image")
    else
        print("\n\nInvalid image")
    end
end

function swapby!(image, booleans, shift)
    dy,dx = shift
    h,w = size(image)
    for y = 1:h, x = 1:w
        if booleans[y,x]
            fy = mod1(y+dy,h)
            fx = mod1(x+dx,w)
            image[y,x], image[fy,fx] = image[fy,fx], image[y,x]
        end
    end
end

# speed(f::Function, args...) = speed(f,1,args...)

function speed(f::Function, times::Int, args...)
    Profile.init(10000000,0.01)
    Profile.clear()
    t=(Float32)[]
    for i in 1:times
        tic()
        @profile f(args...)
        push!(t,toq())
    end
    ProfileView.view()
    println("The function took $(sum(t)/length(t)) seconds on average. ")
end
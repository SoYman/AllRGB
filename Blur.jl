using Images
using OpenCL

const cl = OpenCL

box(image, radius) = imfilter(image, ones(Float32,radius,radius)/(radius^2))
line(image, radius) = imfilter(image, ones(Float32,radius,1)/(radius))

function gaussrgb(image::Image, sigma, iterations)
    image["spatialorder"] = ["y","x"]
    wl, wu, m = solveinteg(sigma, iterations)
    result = convert(Array{RGB{Float32}}, copy(image))
    h,w = size(image)

    for i = 1:iterations
        rad = i-1<m?round(Int, (wl-1)/2):round(Int, (wu-1)/2)
        up = rad+1

        intim = cumsum(cumsum(padarray(result, (up,up),(rad,rad), "circular"),1),2)

        expanded = circshift(intim,(rad,rad)) - circshift(intim,(rad ,-up )) - circshift(intim,(-up, rad)) + circshift(intim,(-up ,-up ))
        
        result = (expanded[up:h+rad,up:w+rad]) / (up+rad)^2
        # view(result)
    end
    return Image(result)
end

function gaussrgb1d(image, sigma, iterations)
    image = Image(image)
    image["spatialorder"] = ["y","x"]
    wl, wu, m = solveinteg(sigma, iterations)
    result = convert(Array{RGB{Float32}}, copy(image))
    h,w = size(image)
    rad,up = 0,0
    for i = 1:iterations*2
        rad = round(Int, i/2)-1<m?round(Int, (wl-1)/2):round(Int, (wu-1)/2)
        up = rad+1

        intim = cumsum(padarray(result, (0,up),(0,rad), "circular"),2)

        expanded = circshift(intim,(0, rad)) - circshift(intim,(0, -up))
        
        result = (expanded[:,up:(i>iterations?h:w)+rad])
        if i == iterations
            result = permutedims(result / ((up+rad)^(iterations)), [2, 1])
        end
    end
    return Image(permutedims(result / ((up+rad)^(iterations)), [2, 1]))
end

function gausslab(image, sigma, iterations)
    image = Image(image)
    image["spatialorder"] = ["y","x"]
    wl, wu, m = solveinteg(sigma, iterations)
    result = convert(Array{Lab{Float32}}, copy(image))
    h,w = size(image)
    rad,up = 0,0
    for i = 1:iterations*2
        rad = round(Int, i/2)-1<m ? round(Int, (wl-1)/2) : round(Int, (wu-1)/2)
        up = rad+1

        intim = cumsum(padarray(result, (0,up),(0,rad), "circular"),2)

        expanded = circshift(intim,(0, rad)) - circshift(intim,(0, -up))
        
        result = (expanded[:,up:(i>iterations?h:w)+rad])
        if i == iterations
            result = transpose(result / ((up+rad)^(iterations)))
        end
    end
    return Image(transpose(result / ((up+rad)^(iterations))))
end

function gaussluv(image, sigma, iterations)
    image = Image(image)
    image["spatialorder"] = ["y","x"]
    wl, wu, m = solveinteg(sigma, iterations)
    result = convert(Array{LUV{Float32}}, copy(image))
    h,w = size(image)
    rad,up = 0,0
    for i = 1:iterations*2
        rad = round(Int, i/2)-1<m?round(Int, (wl-1)/2):round(Int, (wu-1)/2)
        up = rad+1

        intim = cumsum(padarray(result, (0,up),(0,rad), "circular"),2)

        expanded = circshift(intim,(0, rad)) - circshift(intim,(0, -up))
        
        result = (expanded[:,up:(i>iterations?h:w)+rad])
        if i == iterations
            result = transpose(result / ((up+rad)^(iterations)))
        end
    end
    return Image(transpose(result / ((up+rad)^(iterations))))
end

function gausslms(image::Image, sigma, iterations)
    image["spatialorder"] = ["y","x"]
    wl, wu, m = solveinteg(sigma, iterations)
    result = convert(Array{LMS{Float32}}, copy(image))
    h,w = size(image)

    for i = 1:iterations
        rad = i-1<m?round(Int, (wl-1)/2):round(Int, (wu-1)/2)
        up = rad+1

        intim = cumsum(cumsum(padarray(result, (up,up),(rad,rad), "circular"),1),2)

        expanded = circshift(intim,(rad,rad)) - circshift(intim,(rad ,-up )) - circshift(intim,(-up, rad)) + circshift(intim,(-up ,-up ))
        
        result = (expanded[up:h+rad,up:w+rad])
        # view(result)
    end
    return Image(result / (up+rad)^(2*iterations))
end

function integradius(sigma, n)
    wl, wu, m = solveinteg(sigma, n)
    sum([i<m?wl:wu for i = 0:n-1])
end

function solveinteg(sigma, n)
    
    if sigma < 0.8
        prround(Int, "Sigma values below about 0.8 cannot be represented")
        sigma = 0.8
    end
    wIdeal = sqrt(12*sigma^2/n + 1) # Ideal averaging filter width    
    # wl is first odd valued integer less than wIdeal
    wl = round(Int, floor(wIdeal))
    mod(wl,2) == 0 ? wl = wl-1 : 0
    # wu is the next odd value > wl
    wu = wl+2
    # Compute m.  Refer to the tech note for derivation of this formula
    mIdeal = (12*sigma^2 - n*wl^2 - 4*n*wl - 3*n)/(-4*wl - 4)
    m = round(mIdeal)
    
    if m > n || m < 0
        error("calculation of m has failed")
    end
    sigmaActual = sqrt((m*wl^2 + (n-m)*wu^2 - n)/12)
    # prround(Int, "$wl $wu  $m actual sigma $sigmaActual\n")
    return wl, wu, m, sigmaActual
end

# OpenCL stuff

boxblurkernel = "
__kernel void boxblur(__global float *in,
                      __global float *out, 
                      ushort const w,
                      ushort const radius)
{
    const int gid = get_global_id(0) * w;
   
    out[gid] = 0;
    
    for(int i = -radius; i < radius+1; i++){
        out[gid] += in[(i+radius+w)%w+gid];
    }
    out[gid] /=  (radius*2+1);
    for(int i = 1; i < w; i++) {
        out[i+gid] = out[i-1+gid] + (in[(i+radius+1+w)%w+gid] - in[(i-radius+w)%w+gid]) / (radius*2+1);
    }
    return ;
}";

boxblurkernelh = "
__kernel void boxblurh(__global float *in,
                       __global float *out, 
                       ushort const w,
                       ushort const h,
                       ushort const radius)
{
    int gid = get_global_id(0);
   
    out[gid] = 0;
    
    for(int i = -radius; i < radius; i++){
        out[gid] += in[((i+h)%h)*w+gid];
    }
    out[gid] /= (radius*2+1);
    for(int i = 1; i < h; i++) {
        out[i*w+gid] = out[(i-1)*w+gid] + (in[((i+radius+h)%h)*w+gid] - in[((i-radius+h)%h)*w+gid]) / (radius*2+1);
    }
    return ;
}";


function gaussCL(image::Image, sigma, iterations)
    # Calculate convolution radius from sigma
    wl, wu, m = solveinteg(sigma, iterations)

    # Initialise GPU computing context
    ctx   = cl.Context(cl.devices())
    queue = cl.CmdQueue(ctx)

    # Host arrays/buffers
    hostbuffer = convert(Array{Float32}, raw(image))
    out = Array(Float32,size(image))
    h,w = size(image)
    h,w = uint16(h), uint16(w)
    println("host")

    # GPU buffers
    buff1 = cl.Buffer(Float32, ctx, (:rw, :copy), hostbuf=hostbuffer)
    buff2 = cl.Buffer(Float32, ctx, (:rw), length(out))
    println("buffers")

    # Kernel compilation
    program = cl.Program(ctx, source=boxblurkernel) |> cl.build!
    kernel = cl.Kernel(program, "boxblur")

    programh = cl.Program(ctx, source=boxblurkernelh) |> cl.build!
    kernelh = cl.Kernel(programh, "boxblurh")
    println("kernels")

    for i = 1:iterations
        radius = uint16(i-1<m?round(Int, (wl-1)/2):round(Int, (wu-1)/2))
        println("starting pass $i")
        println("$((height(out))), $(width(out)), $radius")
        cl.call(queue, kernel, w, nothing, buff1, buff2, h, radius)
        println("horizontal pass")
        cl.call(queue, kernelh, h, nothing, buff2, buff1, h, w, radius)
        println("vertical pass")
    end
    cl.copy!(queue, out, buff2)

    return (Image(out))
end

# Commented out since it needs an update
# function enqueue_gauss_kernel{T}(queue::cl.CmdQueue, kernel::cl.Kernel, dst::cl.Buffer{T}, src::cl.Buffer{T}, dims)
#     BLOCK_SIZE = 16
#     h,w = dims
#     @assert w % BLOCK_SIZE == 0
#     @assert h % BLOCK_SIZE == 0
#     cl.set_args!(kernel, dst, src, uint32(h), uint32(h))
#     cl.enqueue_kernel(queue, kernel, (h, w), (BLOCK_SIZE, BLOCK_SIZE))

# end

function split_colors(image::Image)
    #function body
end
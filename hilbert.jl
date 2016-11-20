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

@debug function hindexInverse(m, h)
    dim = 3
    e = 0
    d = 0
    p1, p2, p3 = 0,0,0
    # p = [0,0,0]
    # @bp
    w,l,e,d = 0,0,0,0
    for i in reverse(0:m-1)
        w = bit(h, 3*i+2) << 2 | bit(h, 3*i+1) << 1 | bit(h, 3*i)
        l = grayCode(w)
        l = TInverse(e, d, dim, l)
        # for j = 1:dim
        #     p[j] |= bit(l, j-1) << i
        # end
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

function tsbold(i)
    # Split on last 4 bits
    if ((i & 7) == 7)
        # All 4 set
        return 4 + tsb(i >> 4)
    end
    # tsb4 = 
    return (0, 1, 0, 2, 0, 1, 0, 3)[(i & 7)+1]
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


function testHilbertInvertible()
    # @bp
    for i in 0:(2^21-2)
        # print(i," ", hindex(8, hindexInverse(8, i)),"\n")
        if hindex(7, hindexInverse(7, i)) != i
            print("failure on ", i)
            return
        end
    end
end
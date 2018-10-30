# voronoi.jl
# algorithm for creating a voronoi tessellation

using GeometryTypes
using StaticArrays

#include("predicates.jl")

imod(a,b) = mod(a-1,b) + 1

function flip23(teta::Tetrahedron, tetb::Tetrahedron)::NTuple{3, Tetrahedron}
    i::Int = 1
    # Faster to use [teta...] than the raw type??
    a::Int = setdiff([teta...], [tetb...])[1]
    while teta[1] != a
        teta = permute(teta,i)
        i += 1
    end
    e::Int = setdiff([tetb...], [teta...])[1]
    while tetb[4] != e
        tetb = permute(tetb,i)
        i += 1
    end
    (Tetrahedron(teta[2],teta[4],e,a),
     Tetrahedron(teta[4],teta[3],e,a),
     Tetrahedron(teta[3],teta[2],e,a))
 end

function flip23_new(teta::Tetrahedron, tetb::Tetrahedron)::NTuple{3,Tetrahedron}
     a::Int = setdiff([teta...],[tetb...])[1]
     e::Int = setdiff([tetb...],[teta...])[1]
     sharedface::Array{Int,1} = intersect([teta...],[tetb...])
     #sharedface = (findfirst(isequal(a),teta)%2 == 1) ? sharedface : sharedface[end:-1:1]
     (Tetrahedron(sharedface[1],sharedface[3],e,a),
      Tetrahedron(sharedface[3],sharedface[2],e,a),
      Tetrahedron(sharedface[2],sharedface[1],e,a))
  end

function flip32(teta::Tetrahedron, tetb::Tetrahedron, tetc::Tetrahedron)::NTuple{2, Tetrahedron}
    # 3 times faster this way!!
    a::Int, e::Int = intersect([teta...], [tetb...], [tetc...])
    i::Int = 1
    while teta[4] != a
        teta = permute(teta,i)
        i += 1
    end
    while teta[3] != e
        teta = rotate(teta,4)
    end
    b::Int = teta[1]
    d::Int = teta[2]

    while tetb[4] != a
        tetb = permute(tetb,i)
        i += 1
    end
    while tetb[3] != e
        tetb = rotate(tetb,4)
    end

    while tetc[4] != a
        tetc = permute(tetc,i)
        i += 1
    end
    while tetc[3] != e
        tetc = rotate(tetc,4)
    end
    c::Int = (tetb[2] == b) ? tetb[1] : tetc[1]

    (Tetrahedron(a,b,c,d),
     Tetrahedron(b,c,d,e))
 end

function flip32_new(teta::Tetrahedron, tetb::Tetrahedron, tetc::Tetrahedron)::NTuple{2, Tetrahedron}
    a::Int, e::Int = intersect(teta, tetb, tetc)
    sharedface::Array{Int,1} = union(setdiff(teta,[a],[e]),
                                 setdiff(tetb,[a],[e]),
                                 setdiff(tetc,[a],[e]))

    (Tetrahedron(a,sharedface[1],sharedface[2],sharedface[3]),
     Tetrahedron(sharedface[1],sharedface[2],sharedface[3],e))
end

function flip14(tet::Tetrahedron, e::Int)::NTuple{4,Tetrahedron}
    (Tetrahedron(tet[1],tet[2],tet[3],e),
     Tetrahedron(tet[1],tet[3],tet[4],e),
     Tetrahedron(tet[3],tet[2],tet[4],e),
     Tetrahedron(tet[1],tet[4],tet[2],e))
end

function flip41(teta::Tetrahedron, tetb::Tetrahedron, tetc::Tetrahedron, tetd::Tetrahedron)::Tetrahedron
     e::Int = intersect([teta...],[tetb...],[tetc...],[tetd...])[1]
     i:Int = 1
     while teta[4] != e
         teta = permute(teta,i)
         i += 1
     end
     while tetb[4] != e
         tetb = permute(tetb,i)
         i += 1
     end
     d::Int = setdiff([tetb...],[teta...])[1]
     Tetrahedron(teta[1],teta[2],teta[3],d)
 end

# TODO: flip44 

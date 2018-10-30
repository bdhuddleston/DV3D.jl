# voronoi.jl
# algorithm for creating a voronoi tessellation

using GeometryTypes
using StaticArrays

#include("predicates.jl")

imod(a,b) = mod(a-1,b) + 1

function flip23(teta::Tetrahedron, tetb::Tetrahedron)::NTuple{3, Tetrahedron}
    i::Int = 1
    # Faster to use [teta...] than the raw type??
    a::Int = setdiff([teta.vertices...], [tetb.vertices...])[1]
    while teta.vertices[1] != a
        teta = permute(teta,i)
        i += 1
    end
    e::Int = setdiff([tetb.vertices...], [teta.vertices...])[1]
    while tetb.vertices[4] != e
        tetb = permute(tetb,i)
        i += 1
    end
    (Tetrahedron([teta.vertices[2],teta.vertices[4],e,a],
                  teta.neighbors[2],teta.neighbors[4],tetb.neighbors[4],teta.neighbors[1]),
     Tetrahedron([teta.vertices[4],teta.vertices[3],e,a],
                  teta.neighbors[4],teta.neighbors[3],tetb.neighbors[4],teta.neighbors[1]),
     Tetrahedron([teta.vertices[3],teta.vertices[2],e,a],
                  teta.neighbors[3],teta.neighbors[2],tetb.neighbors[4],teta.neighbors[1]))
 end

function flip32(teta::Tetrahedron, tetb::Tetrahedron, tetc::Tetrahedron)::NTuple{2, Tetrahedron}
    # 3 times faster this way!!
    a::Int, e::Int = intersect([teta.vertices...], [tetb.vertices...], [tetc.vertices...])
    i::Int = 1
    while teta.vertices[4] != a
        teta = permute(teta,i)
        i += 1
    end
    while teta.vertices[3] != e
        teta = rotate(teta,4)
    end
    b::Int = teta.vertices[1]
    d::Int = teta.vertices[2]

    while tetb.vertices[4] != a
        tetb = permute(tetb,i)
        i += 1
    end
    while tetb.vertices[3] != e
        tetb = rotate(tetb,4)
    end

    while tetc.vertices[4] != a
        tetc = permute(tetc,i)
        i += 1
    end
    while tetc.vertices[3] != e
        tetc = rotate(tetc,4)
    end
    c::Int = (tetb.vertices[2] == b) ? tetb.vertices[1] : tetc.vertices[1]
    cn::Int = (tetb.vertices[2] == b) ? tetb.neighbors[1] : tetc.neighbors[1]

    (Tetrahedron([a,b,c,d],
                 [teta.neighbors[4],teta.neighbors[1],cn,teta.neighbors[2]]),
     Tetrahedron([b,c,d,e],
                 [teta.neighbors[1],cn,teta.neighbors[2],teta.neighbors[3]]))
 end

function flip14(tet::Tetrahedron, e::Int, offset::Int)::NTuple{4,Tetrahedron}
    (Tetrahedron([tet.vertices[1],tet.vertices[2],tet.vertices[3],e],
                 [offset+3,       offset+2,       offset+4,tet.neighbors[4]]),
     Tetrahedron([tet.vertices[1],tet.vertices[3],tet.vertices[4],e],
                 [offset+3,       offset+4,       offset+1,tet.neighbors[2]]),
     Tetrahedron([tet.vertices[3],tet.vertices[2],tet.vertices[4],e],
                 [offset+4,       offset+2,       offset+1,tet.neighbors[1]]),
     Tetrahedron([tet.vertices[1],tet.vertices[4],tet.vertices[2],e],
                 [offset+3,       offset+1,       offset+2,tet.neighbors[3]]))
end

# TODO: If needed, define a flip41 - not sure that I need the functionality to delete a point
#function flip41(teta::Tetrahedron, tetb::Tetrahedron, tetc::Tetrahedron, tetd::Tetrahedron)::Tetrahedron
#     e::Int = intersect([teta.vertices...],[tetb.vertices...],[tetc.vertices...],[tetd.vertices...])[1]
#     i:Int = 1
#     while teta.vertices[4] != e
#         teta = permute(teta,i)
#         i += 1
#     end
#     while tetb.vertices[4] != e
#         tetb = permute(tetb,i)
#         i += 1
#     end
#     while tetc.vertices[4] != e
#         tetc = permute(tetc,i)
#         i += 1
#     end
#     while tetd.vertices[4] != e
#         tetd = permute(tetd,i)
#         i += 1
#     end
#     d::Int = setdiff([tetb.vertices...],[teta.vertices...])[1]
#     Tetrahedron([teta.vertices[1],teta.vertices[2],teta.vertices[3],d],
#                 [])
#end

function flip44(teta::Tetrahedron, tetb::Tetrahedron, tetc::Tetrahedron, tetd::Tetrahedron)::NTuple{4,Tetrahedron}
    # Tets need to be in the correct order prior to this function!!
    # abcd, aced, fecd, fcbd, where bced is the degenerate plane
    a::Int, b::Int, c::Int, d::Int = teta.vertices
    f::Int = tetc.vertices[1]
    e::Int = tetc.vertices[2]
    (Tetrahedron([a,b,e,d],
                 [teta.neighbors[1],teta.neighbors[2],tetc.neighbors[2],teta.neighbors[4]]),
     Tetrahedron([a,b,c,e],
                 [teta.neighbors[1],teta.neighbors[2],teta.neighbors[3],tetc.neighbors[2]]),
     Tetrahedron([f,d,e,b],
                 [tetc.neighbors[1],teta.neighbors[4],tetc.neighbors[2],teta.neighbors[2]]),
     Tetrahedron([f,e,c,b],
                 [tetc.neighbors[1],tetc.neighbors[2],teta.neighbors[3],teta.neighbors[2]]))
end

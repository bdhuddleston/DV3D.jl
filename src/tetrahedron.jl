# tetrhedron.jl

using GeometryTypes
using StaticArrays


struct Tetrahedron
    vertices::SVector{4,Int}
    neighbors::MVector{4,Int}
end

function permute(tet::Tetrahedron, evenorodd::Int)::Tetrahedron
    if (evenorodd % 2) == 0
        Tetrahedron([tet.vertices[2],tet.vertices[1],tet.vertices[4],tet.vertices[3]],
                          [tet.neighbors[2],tet.neighbors[1],tet.neighbors[4],tet.neighbors[3]])
    else
        Tetrahedron([tet.vertices[4],tet.vertices[3],tet.vertices[2],tet.vertices[1]],
                          [tet.neighbors[4],tet.neighbors[3],tet.neighbors[2],tet.neighbors[1]])
    end
end

function rotate(tet::Tetrahedron, pivot::Int)::Tetrahedron
    if pivot == 1
        Tetrahedron([tet.vertices[1], tet.vertices[4], tet.vertices[2], tet.vertices[3]],
                  [tet.neighbors[1],tet.neighbors[4],tet.neighbors[2],tet.neighbors[3]])
    elseif pivot == 2
        Tetrahedron([tet.vertices[4], tet.vertices[2], tet.vertices[1], tet.vertices[3]],
                  [tet.neighbors[4],tet.neighbors[2],tet.neighbors[1],tet.neighbors[3]])
    elseif pivot == 3
        Tetrahedron([tet.vertices[4], tet.vertices[1], tet.vertices[3], tet.vertices[2]],
                  [tet.neighbors[4],tet.neighbors[1],tet.neighbors[3],tet.neighbors[2]])
    elseif pivot == 4
        Tetrahedron([tet.vertices[3], tet.vertices[1], tet.vertices[2], tet.vertices[4]],
                  [tet.neighbors[3],tet.neighbors[1],tet.neighbors[2],tet.neighbors[4]])
    else
        tet
    end
end

function getface(tet::Tetrahedron, index::Int)::Face{3,Int}
    face::Array{Int,1} = tet.vertices[setdiff(1:4,index)]
    if index % 2 == 1
        Face{3,Int}(face[1],face[3],face[2])
    else
        Face{3,Int}(face...)
    end
end

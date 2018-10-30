# testing

using StaticArrays

abstract type AbstractTet end

struct BoundaryTet <: AbstractTet
end

struct Tetrahedron <: AbstractTet
    vertices::MVector{4,Int}
    neighbors::MVector{4,Int}
end


Tetrahedron() = Tetrahedron([-1,-1,-1,-1],[-1,-1,-1,-1])

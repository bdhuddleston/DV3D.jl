# tetrhedron.jl

# TODO: Define functions to get faces from tets 

const Tetrahedron = SVector{4,Int}

function permute(tet::Tetrahedron, evenorodd::Int)::Tetrahedron
    if (evenorodd % 2) == 0
        tet = Tetrahedron(tet[2],tet[1],tet[4],tet[3])
    else
        tet = Tetrahedron(tet[4],tet[3],tet[2],tet[1])
    end
end

function rotate(tet::Tetrahedron, pivot::Int)::Tetrahedron
    if pivot == 1
        Tetrahedron(tet[1], tet[4], tet[2], tet[3])
    elseif pivot == 2
        Tetrahedron(tet[4], tet[2], tet[1], tet[3])
    elseif pivot == 3
        Tetrahedron(tet[4], tet[1], tet[3], tet[2])
    elseif pivot == 4
        Tetrahedron(tet[3], tet[1], tet[2], tet[4])
    else
        tet
    end
end

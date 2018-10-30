# insertion.jl
# functions for inserting a point into a Delaunay tessellation

#function insertpoint()
# TODO: insertpoint! (overload push!() ?)
# INPUT: DT object, point
# OUTPUT: nothing
#
# Walk to containing tet
# Add point to DT points
# pop containing tet
# call a flip14 on the containing tet and append tets to DT tets
# loop through old tet neighbors, updating neighbor information (subroutine)
#
# Push new tets on stack
# while stack is not empty
#   pop a tet from the stack -> pabc
#   Get the neighbor opposite the new point -> abcd
#   if insphere of tet pabc and d is positive
#       Flip the tets (subroutine)
#
#


# TODO: update neighbors
# INPUT: old tet neighbors, old tet index, new tets
# OUTPUT: nothing

# TODO: flip neighboring tets 

# walk.jl
# implement walk method to locate a point

# TODO: Walk function
# INPUT: index of starting tet, point to locate, Array of points, Array of tets
# OUTPUT: index of tet containing point
#
# While we haven't found the containing tet (new index does not match last index)
#   update the current tet with STEP or the current index
#

# TODO: Step function
# INPUT: Current tet, points of current tet, point to locate
# OUTPUT: Index of next tet
#
# Iterate through the faces
#   if orient of the face and the target point is negative
#       return the neighbor index associated with that face
# if it falls through, return nothing
#

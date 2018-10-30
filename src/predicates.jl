# fast robust geometric predicates
# Converted into Julia from C
# Original functions written by
#  Jonathan Richard Shewchuk
#  School of Computer Science
#  Carnegie Mellon University
#  5000 Forbes Avenue
#  Pittsburgh, Pennsylvania  15213-3891
#  jrs@cs.cmu.edu

using GeometryTypes


"""
original comment:
/* Which of the following two methods of finding the absolute values is      */
/*   fastest is compiler-dependent.  A few compilers can inline and optimize */
/*   the fabs() call; but most will incur the overhead of a function call,   */
/*   which is disastrously slow.  A faster way on IEEE machines might be to  */
/*   mask the appropriate bit, but that's difficult to do in C.              */
This function definition attempts to do the same thing as a C-Define
"""
@inline absolute(a) = ((a) >= 0.0 ? (a) : -(a))
# TODO: should also try the built in absolute
# @inline absolute(a) = abs(a)

# Skipping tails, condensing to one function (?)
#@inline function fast_two_sum_tail(a::Float64, b::Float64, x::Float64)::Tuple{Float64,Float64}
#    bvirt::Float64 = x - a
#    x, b - bvirt
#end

@inline function fast_two_sum(a::Float64, b::Float64)::Tuple{Float64,Float64}
    x::Float64 = a + b
    bvirt::Float64 = x - a
    x, b - bvirt
end

#@inline function fast_two_diff_tail(a::Float64, b::Float64, x::Float64)

const resulterrbound = Ref{Float64}()
const ccwerrboundA = Ref{Float64}()
const ccwerrboundB = Ref{Float64}()
const ccwerrboundC = Ref{Float64}()
const o3derrboundA = Ref{Float64}()
const o3derrboundB = Ref{Float64}()
const o3derrboundC = Ref{Float64}()
const iccerrboundA = Ref{Float64}()
const iccerrboundB = Ref{Float64}()
const iccerrboundC = Ref{Float64}()
const isperrboundA = Ref{Float64}()
const isperrboundB = Ref{Float64}()
const isperrboundC = Ref{Float64}()
"""
/*****************************************************************************/
/*                                                                           */
/*  exactinit()   Initialize the variables used for exact arithmetic.        */
/*                                                                           */
/*  `epsilon' is the largest power of two such that 1.0 + epsilon = 1.0 in   */
/*  floating-point arithmetic.  `epsilon' bounds the relative roundoff       */
/*  error.  It is used for floating-point error analysis.                    */
/*                                                                           */
/*  `splitter' is used to split floating-point numbers into two half-        */
/*  length significands for exact multiplication.                            */
/*                                                                           */
/*  I imagine that a highly optimizing compiler might be too smart for its   */
/*  own good, and somehow cause this routine to fail, if it pretends that    */
/*  floating-point arithmetic is too much like real arithmetic.              */
/*                                                                           */
/*  Don't change this routine unless you fully understand it.                */
/*                                                                           */
/*****************************************************************************/
"""
function exactinit()::Nothing
    every_other::Bool = false
    epsilon::Float64 = 1.0
    splitter::Float64 = 1.0
    check::Float64 = 1.0
    #/* Repeatedly divide `epsilon' by two until it is too small to add to    */
    #/*   one without causing roundoff.  (Also check if the sum is equal to   */
    #/*   the previous sum, for machines that round up instead of using exact */
    #/*   rounding.  Not that this library will work on such machines anyway. */
    lastcheck = check;
    epsilon *= 0.5
    splitter *= 2.0
    check = 1.0 + epsilon
    while ((check != 1.0) && (check != lastcheck))
        lastcheck = check;
        epsilon *= 0.5
        if (every_other)
            splitter *= 2.0
        end
        every_other = !every_other
        check = 1.0 + epsilon
    end
    splitter += 1.0

    #/* Error bounds for orientation and incircle tests. */
    global resulterrbound[] = (3.0 + 8.0 * epsilon) * epsilon
    global ccwerrboundA[] = (3.0 + 16.0 * epsilon) * epsilon
    global ccwerrboundB[] = (2.0 + 12.0 * epsilon) * epsilon
    global ccwerrboundC[] = (9.0 + 64.0 * epsilon) * epsilon * epsilon
    global o3derrboundA[] = (7.0 + 56.0 * epsilon) * epsilon
    global o3derrboundB[] = (3.0 + 28.0 * epsilon) * epsilon
    global o3derrboundC[] = (26.0 + 288.0 * epsilon) * epsilon * epsilon
    global iccerrboundA[] = (10.0 + 96.0 * epsilon) * epsilon
    global iccerrboundB[] = (4.0 + 48.0 * epsilon) * epsilon
    global iccerrboundC[] = (44.0 + 576.0 * epsilon) * epsilon * epsilon
    global isperrboundA[] = (16.0 + 224.0 * epsilon) * epsilon
    global isperrboundB[] = (5.0 + 72.0 * epsilon) * epsilon
    global isperrboundC[] = (71.0 + 1408.0 * epsilon) * epsilon * epsilon
    nothing
end


"""
/*               Return a positive value if the point pd lies below the      */
/*               plane passing through pa, pb, and pc; "below" is defined so */
/*               that pa, pb, and pc appear in counterclockwise order when   */
/*               viewed from above the plane.  Returns a negative value if   */
/*               pd lies above the plane.  Returns zero if the points are    */
/*               coplanar.  The result is also a rough approximation of six  */
/*               times the signed volume of the tetrahedron defined by the   */
/*               four points.                                                */
"""
function orient3d(pa::Point3f0, pb::Point3f0, pc::Point3f0, pd::Point3f0)::Float64
    ad::Point3f0 = pa .- pd
    bd::Point3f0 = pb .- pd
    cd::Point3f0 = pc .- pd

    bdxcdy::Float64 = bd[1] * cd[2]
    cdxbdy::Float64 = cd[1] * bd[2]
    cdxady::Float64 = cd[1] * ad[2]
    adxcdy::Float64 = ad[1] * cd[2]
    adxbdy::Float64 = ad[1] * bd[2]
    bdxady::Float64 = bd[1] * ad[2]

    det::Float64 = adz * (bdxcdy - cdxbdy) +
                   bdz * (cdxady - adxcdy) +
                   cdz * (adxbdy - bdxady)
    permanent = (absolute(bdxcdy) + absolute(cdxbdy)) * absolute(adz) +
                (absolute(cdxady) + absolute(adxcdy)) * absolute(bdz) +
                (absolute(adxbdy) + absolute(bdxady)) * absolute(cdz)
    # TODO: define predetermined error bounds
    errbound = o3derrboundA * permanent
    if (det > errbound) || (-det > errbound)
        return det
    else
        println("WARNING: possible errors in floating point arithmetic!!!")
        return det
    end
    # TODO: implement robust predicates
    #return orient3dadapt(pa, pb, pc, pd, permanent)
end


"""
/*               Return a positive value if the point pe lies inside the     */
/*               sphere passing through pa, pb, pc, and pd; a negative value */
/*               if it lies outside; and zero if the five points are         */
/*               cospherical.  The points pa, pb, pc, and pd must be ordered */
/*               so that they have a positive orientation (as defined by     */
/*               orient3d()), or the sign of the result will be reversed.    */
"""
function insphere(pa::Point3f0, pb::Point3f0, pc::Point3f0, pd::Point3f0, pe::Point3f0)::Float64
    ae::Point3f0 = pa .- pe
    be::Point3f0 = pb .- pe
    ce::Point3f0 = pc .- pe
    de::Point3f0 = pd .- pe

    aexbey = ae[1] * be[2]
    bexaey = be[1] * ae[2]
    ab = aexbey - bexaey
    bexcey = be[1] * ce[2]
    cexbey = ce[1] * be[2]
    bc = bexcey - cexbey
    cexdey = ce[1] * de[2]
    dexcey = de[1] * ce[2]
    cd = cexdey - dexcey
    dexaey = de[1] * ae[2]
    aexdey = ae[1] * de[2]
    da = dexaey - aexdey

    aexcey = ae[1] * ce[2]
    cexaey = ce[1] * ae[2]
    ac = aexcey - cexaey
    bexdey = be[1] * de[2]
    dexbey = de[1] * be[2]
    bd = bexdey - dexbey

    abc = ae[3] * bc - be[3] * ac + ce[3] * ab
    bcd = be[3] * cd - ce[3] * bd + de[3] * bc
    cda = ce[3] * da + de[3] * ac + ae[3] * cd
    dab = de[3] * ab + ae[3] * bd + be[3] * da

    alift = sum(ae .* ae)
    blift = sum(be .* be)
    clift = sum(ce .* ce)
    dlift = sum(de .* de)

    det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd)
    permanent = ((absolute(cexdey) + absolute(dexcey)) * absolute(bez) +
                 (absolute(dexbey) + absolute(bexdey)) * absolute(cez) +
                 (absolute(bexcey) + absolute(cexbey)) * absolute(dez)) * alift +
                ((absolute(dexaey) + absolute(aexdey)) * absolute(cez) +
                 (absolute(aexcey) + absolute(cexaey)) * absolute(dez) +
                 (absolute(cexdey) + absolute(dexcey)) * absolute(aez)) * blift +
                ((absolute(aexbey) + absolute(bexaey)) * absolute(dez) +
                 (absolute(bexdey) + absolute(dexbey)) * absolute(aez) +
                 (absolute(dexaey) + absolute(aexdey)) * absolute(bez)) * clift +
                ((absolute(bexcey) + absolute(cexbey)) * absolute(aez) +
                 (absolute(cexaey) + absolute(aexcey)) * absolute(bez) +
                 (absolute(aexbey) + absolute(bexaey)) * absolute(cez)) * dlift
    errbound = isperrboundA * permanent;
    if (det > errbound) || (-det > errbound)
        return det
    else
        println("WARNING: possible error in floating point arithmetic!!!")
        return det
    end
    # TODO: implement adaptive predicates
end

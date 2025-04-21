using KiteUtils
using KiteModels


function print_inertia_matrix(Ixx, Ixy, Ixz, Iyy, Iyz, Izz)
    println("Inertia matrix [kgm²]:")
    println(" Ixx Ixy Ixz: [$Ixx $Ixy $Ixz] ")
    println(" Ixy Iyy Iyz: [$Ixy $Iyy $Iyz] ")
    println(" Ixz Iyz Izz: [$Ixz $Iyz $Izz] ")
end


function print_settings(include_kcu::Bool=true, around_kcu::Bool=false)
    out = ""

    if include_kcu
        out *= "Calculating the inertia for kcu and kite"
    else
        out *= "Calculating the inertia for kite alone"
    end

    if around_kcu
        out *= " around the kcu."
    else
        out *= " around the center of mass."
    end 

    println(out)
end


SETFILE = "system.yaml"
INCL_KCU = false
AROUND_KCU = false

SET = deepcopy(load_settings(SETFILE))
TEST_KCU = KCU(SET)
S = KPS4(TEST_KCU)

print_settings(INCL_KCU, AROUND_KCU)
ixx, ixy, ixz, iyy, iyz, izz = calculate_rotational_inertia!(S, INCL_KCU, AROUND_KCU)
print_inertia_matrix(ixx, ixy, ixz, iyy, iyz, izz)
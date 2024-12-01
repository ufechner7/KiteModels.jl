using KiteUtils

function calculate_inertia_for_setting(settings_file::String, include_kcu::Bool=true, around_kcu::Bool=false)
    set_data_path("data")
    set = deepcopy(load_settings(settings_file))

    points = KiteUtils.get_particles(set.height_k, set.h_bridle, set.width, set.m_k, [0, 0, 0], [0, 0, -1], [10, 0, 0])
    
    pos_matrix = [points[begin+1] points[begin+2] points[begin+3] points[begin+4] points[begin+5]]
    X = pos_matrix[begin, :]
    Y = pos_matrix[begin+1, :]
    Z = pos_matrix[begin+2, :]

    k2 = set.rel_top_mass * (1.0 - set.rel_nose_mass)
    k3 = 0.5 * (1.0 - set.rel_top_mass) * (1.0 - set.rel_nose_mass)
    M = [set.kcu_mass, set.rel_nose_mass * set.mass, k2 * set.mass, k3 * set.mass, k3 * set.mass]
    
    if !include_kcu
        X = X[begin+1:end]
        Y = Y[begin+1:end]
        Z = Z[begin+1:end]
        M = M[begin+1:end]
    end

    if around_kcu
        Ixx, Ixy, Ixz, Iyy, Iyz, Izz = calculate_rotational_inertia(X, Y, Z, M, false, points[begin+1])
    else
        Ixx, Ixy, Ixz, Iyy, Iyz, Izz = calculate_rotational_inertia(X, Y, Z, M)
    end

    Ixx, Ixy, Ixz, Iyy, Iyz, Izz
end


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


SETFILE = "system_v9.yaml"
INCL_KCU = false
ARROUND_KCU = false

print_settings(INCL_KCU, ARROUND_KCU)
ixx, ixy, ixz, iyy, iyz, izz = calculate_inertia_for_setting(SETFILE, INCL_KCU, ARROUND_KCU)
print_inertia_matrix(ixx, ixy, ixz, iyy, iyz, izz)

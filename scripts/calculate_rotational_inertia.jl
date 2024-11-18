using KiteUtils

function calculate_rotational_inertia(X::Vector, Y::Vector, Z::Vector, M::Vector, around_center_of_mass::Bool=true, 
    rotation_point::Vector=[0, 0, 0])
    @assert size(X) == size(Y) == size(Z) == size(M)
    
    if around_center_of_mass
        # First loop to determine the center of mass
        x_com = y_com = z_com = m_total = 0.0
        for (x, y, z, m) in zip(X, Y, Z, M)
            x_com += x * m
            y_com += y * m
            z_com += z * m
            m_total += m 
        end

        x_com = x_com / m_total
        y_com = y_com / m_total
        z_com = z_com / m_total
    else
        x_com = rotation_point[begin]
        y_com = rotation_point[begin+1]
        z_com = rotation_point[begin+2]
    end

    Ixx = Ixy = Ixz = Iyy = Iyz = Izz = 0

    # Second loop using the distance between the point and the center of mass
    for (x, y, z, m) in zip(X .- x_com, Y .- y_com, Z .- z_com, M)
        Ixx += m * (y^2 + z^2)
        Iyy += m * (x^2 + z^2)
        Izz += m * (x^2 + y^2)

        Ixy += m * x * y
        Ixz += m * x * z
        Iyz += m * y * z
    end
    
    Ixx, Ixy, Ixz, Iyy, Iyz, Izz
end


function calculate_intertia_for_setting(settings_file::String, include_kcu::Bool=true, around_kcu::Bool=false)
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
INCLKCU = false
ARROUNDKCU = false

print_settings(INCLKCU, ARROUNDKCU)
IXX, IXY, IXZ, IYY, IYZ, IZZ = calculate_intertia_for_setting(SETFILE, INCLKCU, ARROUNDKCU)
print_inertia_matrix(IXX, IXY, IXZ, IYY, IYZ, IZZ)

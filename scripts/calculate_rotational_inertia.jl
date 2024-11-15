function calculate_rotational_inertia(X::Vector, Y::Vector, Z::Vector, M::Vector; around_center_of_mass=true, rotation_point=[0, 0, 0])
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
    
    [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz]
end
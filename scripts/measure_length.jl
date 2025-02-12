using Interpolations
using Plots
using Serialization

# Read vertices from the .obj file
function read_vertices(filename)
    vertices = []
    open(filename) do file
        for line in eachline(file)
            if startswith(line, "v ") && !startswith(line, "vt") && !startswith(line, "vn")
                parts = split(line)
                x = parse(Float64, parts[2])
                y = parse(Float64, parts[3])
                z = parse(Float64, parts[4])
                push!(vertices, (x, y, z))
            end
        end
    end
    return vertices
end

# Create interpolations for max and min x coordinates
function create_interpolations(vertices; gamma_range=-π:0.01:π, box_center_z=-2.45)
    vz_centered = [v[3] - box_center_z for v in vertices]
    
    max_xs = Float64[]
    min_xs = Float64[]
    gammas = Float64[]
    
    for gamma in gamma_range
        cos_g = cos(gamma)
        sin_g = sin(gamma)
        x_inside = Float64[]
        
        for (i, v) in enumerate(vertices)
            # Rotate y coordinate to check box containment
            rotated_y = v[2] * cos_g + vz_centered[i] * sin_g
            if abs(rotated_y) ≤ 0.1
                push!(x_inside, v[1])
            end
        end
        
        if !isempty(x_inside)
            push!(gammas, gamma)
            push!(max_xs, maximum(x_inside))
            push!(min_xs, minimum(x_inside))
        end
    end
    
    # Create sorted interpolations
    order = sortperm(gammas)
    sorted_gammas = gammas[order]
    sorted_max = max_xs[order]
    sorted_min = min_xs[order]
    
    itp_max = linear_interpolation(sorted_gammas, sorted_max)
    itp_min = linear_interpolation(sorted_gammas, sorted_min)
    
    return (itp_max, itp_min, sorted_gammas, sorted_max, sorted_min)
end

# Main execution
vertices = read_vertices("kite.obj")
(itp_max, itp_min, gammas, max_xs, min_xs) = create_interpolations(vertices)
serialize(se.polar_file, (itp_max, itp_min, gammas, max_xs, min_xs))
nothing
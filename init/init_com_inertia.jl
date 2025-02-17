# Read faces from the .obj file
function read_faces(filename)
    vertices = []
    faces = []
    
    open(filename) do file
        for line in eachline(file)
            if startswith(line, "v ") && !startswith(line, "vt") && !startswith(line, "vn")
                parts = split(line)
                x = parse(Float64, parts[2])
                y = parse(Float64, parts[3])
                z = parse(Float64, parts[4])
                push!(vertices, [x, y, z])
            elseif startswith(line, "f ")
                parts = split(line)
                # Handle both f v1 v2 v3 and f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3 formats
                indices = map(p -> parse(Int, split(p, '/')[1]), parts[2:4])
                push!(faces, indices)
            end
        end
    end
    return vertices, faces
end

# Calculate center of mass for a triangular surface mesh
function calculate_com(vertices, faces)
    area_total = 0.0
    com = zeros(3)
    
    for face in faces
        v1 = vertices[face[1]]
        v2 = vertices[face[2]]
        v3 = vertices[face[3]]
        
        # Calculate triangle area and centroid
        normal = cross(v2 - v1, v3 - v1)
        area = norm(normal) / 2
        centroid = (v1 + v2 + v3) / 3
        
        area_total += area
        com += area * centroid
    end
    
    return com / area_total
end

# Calculate inertia tensor for a triangular surface mesh with given mass
function calculate_inertia_tensor(vertices, faces, mass, com)
    # Initialize inertia tensor
    I = zeros(3, 3)
    total_area = 0.0
    
    for face in faces
        v1 = vertices[face[1]] .- com
        v2 = vertices[face[2]] .- com
        v3 = vertices[face[3]] .- com
        
        # Calculate triangle area
        normal = cross(v2 - v1, v3 - v1)
        area = norm(normal) / 2
        total_area += area
        
        # Calculate contribution to inertia tensor
        for i in 1:3
            for j in 1:3
                # Vertices relative to center of mass
                points = [v1, v2, v3]
                
                # Calculate contribution to inertia tensor
                for p in points
                    if i == j
                        # Diagonal terms
                        I[i,i] += area * (sum(p.^2) - p[i]^2)
                    else
                        # Off-diagonal terms
                        I[i,j] -= area * (p[i] * p[j])
                    end
                end
            end
        end
    end
    
    # Scale by mass/total_area to get actual inertia tensor
    return (mass / total_area) * I / 3
end

function calculate_kite_properties(filename, mass)
    vertices, faces = read_faces(filename)
    com = calculate_com(vertices, faces)
    inertia_tensor = calculate_inertia_tensor(vertices, faces, mass, com)
    
    println("Center of Mass:")
    println(com)
    println("\nInertia Tensor:")
    display(inertia_tensor)
    
    return com, inertia_tensor
end


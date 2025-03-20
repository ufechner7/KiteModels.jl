# Test the euler formulas with different frames
using Test
using LinearAlgebra
using KiteModels

@testset "calc_inertia2 tests" begin
    @testset "Standard diagonal inertia tensor" begin
        # Create a diagonal inertia tensor
        I_b_tensor = Diagonal([10.0, 20.0, 30.0])
        
        # Call the function
        I_p, I_b, R_b_p, Q_p_b = KiteModels.calc_inertia2(I_b_tensor)
        
        # For a diagonal tensor, the principal frame should already be aligned with body frame
        @test I_p ≈ [10.0, 20.0, 30.0]
        @test I_b ≈ [10.0, 20.0, 30.0]
        @test R_b_p ≈ I(3)  # Identity matrix
        @test Q_p_b ≈ [1.0, 0.0, 0.0, 0.0]  # Identity quaternion
    end
    
    @testset "Non-diagonal inertia tensor" begin
        # Create non-diagonal inertia tensor
        # This represents a 45° rotation around the y-axis
        theta = π/4
        R = [cos(theta) 0 sin(theta); 
             0 1 0; 
            -sin(theta) 0 cos(theta)]
        
        I_diag = Diagonal([10.0, 20.0, 30.0])
        I_rotated = R * I_diag * R'
        
        # Call the function
        I_p, I_b, R_b_p, Q_p_b = KiteModels.calc_inertia2(I_rotated)
        
        # Check diagonal elements of original tensor
        @test I_b ≈ [I_rotated[1,1], I_rotated[2,2], I_rotated[3,3]]
        
        # Test that principal moments are sorted
        @test I_p[1] ≤ I_p[2] ≤ I_p[3]
        
        # Test that eigenvalues are correct (might be reordered)
        @test sort(I_p) ≈ sort([10.0, 20.0, 30.0])
        
        # Check that R_b_p properly diagonalizes the tensor
        I_p_test = R_b_p' * I_rotated * R_b_p
        @test isapprox(I_p_test, Diagonal(I_p), atol=1e-10)
        
        # Check that R_b_p is a proper rotation matrix
        @test isapprox(det(R_b_p), 1.0, atol=1e-10)
        @test isapprox(R_b_p * R_b_p', I(3), atol=1e-10)
        
        # Check quaternion matches rotation matrix
        R_from_q = KiteModels.quaternion_to_rotation_matrix(Q_p_b)'
        @test isapprox(R_from_q, R_b_p, atol=1e-10)
    end
    
    @testset "Alignment with standard basis" begin
        # Create a rotation around a non-standard axis that results in eigenvectors
        # that aren't aligned with standard basis
        axis = normalize([1.0, 2.0, 3.0])
        theta = π/3
        
        # Rodrigues' rotation formula for arbitrary axis
        K = [0 -axis[3] axis[2]; axis[3] 0 -axis[1]; -axis[2] axis[1] 0]
        R = I(3) + sin(theta) * K + (1 - cos(theta)) * (K^2)
        
        I_diag = Diagonal([5.0, 10.0, 15.0])
        I_rotated = R * I_diag * R'
        
        # Call the function
        I_p, I_b, R_b_p, Q_p_b = KiteModels.calc_inertia2(I_rotated)
        
        # Test that the function preserves proper diagonalization
        I_p_test = R_b_p' * I_rotated * R_b_p
        @test isapprox(I_p_test, Diagonal(I_p), atol=1e-10)
        
        # Test that eigenvalues match
        @test sort(I_p) ≈ sort(diag(I_diag))
        
        # Check for proper rotation matrix properties
        @test isapprox(det(R_b_p), 1.0, atol=1e-10)
        @test isapprox(R_b_p * R_b_p', I(3), atol=1e-10)
        
        # Check that eigenvectors are maximally aligned with standard basis
        # We can test this by checking if any other permutation would give better alignment
        
        # Function to calculate alignment score for a given permutation
        function alignment_score(evecs, perm)
            score = 0.0
            for i in 1:3
                score += abs(evecs[:,i] ⋅ [j==perm[i] ? 1.0 : 0.0 for j in 1:3])
            end
            return score
        end
        
        # Get original eigenvectors
        eigen_result = eigen(I_rotated)
        p = sortperm(eigen_result.values)
        evecs = eigen_result.vectors[:, p]
        
        # Test all possible permutations
        permutations = [
            [1,2,3], [1,3,2], [2,1,3], 
            [2,3,1], [3,1,2], [3,2,1]
        ]
        
        our_alignment = sum(diag(abs.(R_b_p' * I(3))))
        best_alignment = maximum([alignment_score(evecs, perm) for perm in permutations])
        
        # Our alignment should be at least as good as the best permutation
        @test our_alignment ≥ best_alignment - 1e-10
    end
    
    @testset "Right-handedness preservation" begin
        # Create a left-handed inertia tensor
        # First define a proper rotation
        R = [0 1 0; 0 0 1; 1 0 0]  # Cyclic permutation
        I_diag = Diagonal([10.0, 20.0, 30.0])
        I_rotated = R * I_diag * R'
        
        # Make it left-handed by flipping a column
        R_left = copy(R)
        R_left[:,3] *= -1
        I_left = R_left * I_diag * R_left'
        
        # Call the function
        I_p, I_b, R_b_p, Q_p_b = KiteModels.calc_inertia2(I_left)
        
        # Check that we get a right-handed coordinate system
        @test det(R_b_p) ≈ 1.0
        
        # Check that the tensor is properly diagonalized
        I_p_test = R_b_p' * I_left * R_b_p
        @test isapprox(I_p_test, Diagonal(I_p), atol=1e-10)
    end
    
    @testset "Complex inertia tensor" begin
        # Create a complex non-symmetric inertia tensor
        I_complex = [
            10.0 -2.0 -3.0;
            -2.0 15.0 -1.0;
            -3.0 -1.0 20.0
        ]
        
        # Make it symmetric (physically valid inertia tensor must be symmetric)
        I_complex = 0.5 * (I_complex + I_complex')
        
        # Call the function
        I_p, I_b, R_b_p, Q_p_b = KiteModels.calc_inertia2(I_complex)
        
        # Check that diagonal elements are extracted correctly
        @test I_b ≈ [I_complex[1,1], I_complex[2,2], I_complex[3,3]]
        
        # Test that the function diagonalizes the tensor
        I_p_test = R_b_p' * I_complex * R_b_p
        @test isapprox(I_p_test, Diagonal(I_p), atol=1e-10)
        
        # The diagonal elements should match eigenvalues
        eigen_vals = eigvals(I_complex)
        @test sort(I_p) ≈ sort(eigen_vals)
        
        # Verify that eigenvalues are in ascending order for principal moments
        @test I_p[1] ≤ I_p[2] ≤ I_p[3]
        
        # Check rotation matrix properties
        @test isapprox(det(R_b_p), 1.0, atol=1e-10)
        @test isapprox(R_b_p * R_b_p', I(3), atol=1e-10)
    end
end

@testset "Frame transformation and moment calculation tests" begin
    @testset "Body to principal frame moment transformations" begin
        # Create a rotation matrix from body to principal frame (45° around y axis)
        theta = π/4
        R_b_p = [cos(theta) 0   sin(theta); 
                 0          1   0; 
                -sin(theta) 0   cos(theta)]
        
        # Define inertia tensor in body frame (diagonal)
        I_b_diag = [10.0, 20.0, 30.0]
        I_b_tensor = Diagonal(I_b_diag)
        
        # Transform inertia tensor to principal frame
        I_p_tensor = R_b_p * I_b_tensor * R_b_p'
        
        # The inertia tensor in principal frame is no longer diagonal!
        
        # Define angular velocity (zero for this test)
        ω_p = zeros(3)
        
        # Apply an external moment in body frame
        moment_b = [5.0, -3.0, 2.0]
        
        # Transform moment to principal frame
        moment_p = R_b_p * moment_b
        
        # Calculate angular acceleration in principal frame
        # Note: We need to use the full tensor equation here since I_p is not diagonal
        α_p = I_p_tensor \ (moment_p + cross(ω_p, I_p_tensor * ω_p))
        
        # Transform acceleration back to body frame
        α_b = R_b_p' * α_p
        
        # The relationship should now hold
        expected_α_b = I_b_tensor \ moment_b  # Equivalent to moment_b ./ I_b_diag
        
        # Validate results
        @test α_b ≈ expected_α_b
    end
    
    # Test quaternion and rotation matrix consistency
    @testset "Quaternion and rotation matrix consistency" begin
        # Create a rotation matrix representing a 30° rotation around the x-axis
        theta = π/6
        R = [1 0 0; 
             0 cos(theta) -sin(theta); 
             0 sin(theta) cos(theta)]
        
        # Convert to quaternion
        q = KiteModels.rotation_matrix_to_quaternion(R)
        
        # Convert back to rotation matrix
        R2 = KiteModels.quaternion_to_rotation_matrix(q)
        
        # Verify the round-trip conversion is accurate
        @test R ≈ R2 atol=1e-10
        
        # Test quaternion multiplication and rotation
        # Create two rotations: 90° around x, then 90° around y
        # The equivalent rotation matrices
        R1 = [1 0 0; 0 0 -1; 0 1 0]
        R2 = [0 0 1; 0 1 0; -1 0 0]
        q1 = KiteModels.rotation_matrix_to_quaternion(R1)
        q2 = KiteModels.rotation_matrix_to_quaternion(R2)
        
        # Multiply quaternions
        q_combined = KiteModels.quaternion_multiply(q2, q1)
        
        # Combined rotation matrix (apply R1 then R2)
        R_combined = R2 * R1
        
        # Convert combined quaternion to rotation matrix
        R_from_q = KiteModels.quaternion_to_rotation_matrix(q_combined)
        
        # Verify the quaternion multiplication gives the same result as matrix multiplication
        @test R_from_q ≈ R_combined atol=1e-10
    end
end
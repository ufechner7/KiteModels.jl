using Test, LinearAlgebra, KiteUtils, VortexStepMethod
using KiteModels
using Statistics

old_path = get_data_path()
package_data_path = joinpath(dirname(dirname(pathof(KiteModels))), "data")
temp_data_path = joinpath(tempdir(), "data")
Base.Filesystem.cptree(package_data_path, temp_data_path; force=true)
set_data_path(temp_data_path)

# Testing tolerance
const TOL = 1e-5
const BUILD_SYS = true

@testset verbose = true "RamAirKite MTK Model Tests" begin
    # Initialize model
    set = se("system_ram.yaml")
    set.segments = 3
    set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
    set.quasi_static = false
    set.physical_model = "ram"

    @info "Creating s:"
    @time s = RamAirKite(set)

    measure = Measurement()
    s.set.abs_tol = 1e-2
    s.set.rel_tol = 1e-2

    # Initialize at elevation
    measure.sphere_pos .= deg2rad.([80.0 80.0; 1.0 -1.0])
    
    @testset "Model Initialization Chain" begin
        if BUILD_SYS
            # Delete existing problem file to force init!
            @info "Data path: $(get_data_path())"
            prob_path = joinpath(get_data_path(), KiteModels.get_prob_name(s.set))
            if isfile(prob_path)
                @info "Removing existing serialized problem from $prob_path to test full initialization"
                rm(prob_path)
            end
            
            # 1. First time initialization - should create new model
            @info "Testing initial init! (should create new model)..."
            @time KiteModels.init_sim!(s, measure; prn=true)
            
            # Check that serialization worked
            @test isfile(prob_path)
            
            # Check initialization results
            @test !isnothing(s.integrator)
            @test !isnothing(s.sys)
            @test !isnothing(s.point_system)
        end
        s.integrator = nothing
        s.sys = nothing
        
        # Keep references to first integrator and point system
        first_integrator_ptr = objectid(s.integrator)
        first_point_system_ptr = objectid(s.point_system)
            
        # 2. First reinit! - should load from serialized file
        @info "Testing first reinit! (should load serialized file)..."
        @time KiteModels.reinit!(s, measure; prn=true, reload=false)
        next_step!(s)
        
        # Check that it's a new integrator
        second_integrator_ptr = objectid(s.integrator)
        second_point_system_ptr = objectid(s.point_system)
        @test first_integrator_ptr != second_integrator_ptr
        @test first_point_system_ptr == second_point_system_ptr
            
        # 3. Second reinit! - should reuse existing integrator
        @info "Testing second reinit! (should reuse integrator)..."
        @time KiteModels.reinit!(s, measure; prn=true, reload=false)
        
        # This should create a new point system but reuse the existing integrator
        third_integrator_ptr = objectid(s.integrator)
        third_point_system_ptr = objectid(s.point_system)
        @test second_integrator_ptr == third_integrator_ptr # Should be the same 
        @test second_point_system_ptr == third_point_system_ptr
            
        # Get positions from various sources
        pos_integrator, _, _, _, _, _, _, _, _, _, _ = s.get_state(s.integrator)
        sys_state = KiteModels.SysState(s)
        
        # Check dimension consistency
        @test size(pos_integrator, 2) == length(s.point_system.points)
        @test length(sys_state.X) == length(s.point_system.points)

        # Compare positions in different representations
        for (i, point) in enumerate(s.point_system.points)
            # Points' world positions should match integrator positions
            point_pos = point.pos_w
            integ_pos = pos_integrator[:, i]
            sys_state_pos = [sys_state.X[i], sys_state.Y[i], sys_state.Z[i]]
            
            @test isapprox(norm(point_pos), norm(integ_pos), rtol=1e-2)
            @test isapprox(norm(sys_state_pos), norm(integ_pos), rtol=1e-2)
            
            # Positions should not be zero (except ground points)
            if point.type != KiteModels.WINCH  # Skip ground points which might be at origin
                @test norm(point_pos) > 0.1
                @test norm(integ_pos) > 0.1
                @test norm(sys_state_pos) > 0.1
            end
        end
    end
    
    @testset "State Consistency" begin
        KiteModels.reinit!(s, measure, prn=true, reload=false)
    
        # Check quaternion normalization
        _, _, Q_b_w, elevation, azimuth, _, _, _, _, _, _ = s.get_state(s.integrator)
        @test isapprox(norm(Q_b_w), 1.0, atol=TOL)
    
        # Check elevation matches measurement
        @test isapprox(elevation, measure.elevation, atol=1e-2)
    
        # Change measurement and reinitialize
        old_elevation = measure.elevation
        measure.sphere_pos[1,:] .= deg2rad(85.0)
        KiteModels.reinit!(s, measure, prn=true, reload=false)
    
        # Get new state
        _, _, _, elevation_new, _, _, _, _, _, _, _ = s.get_state(s.integrator)
    
        # Verify state changed according to measurement
        @test !isapprox(elevation_new, old_elevation, atol=1e-2)
        @test isapprox(elevation_new, deg2rad(85.0), atol=1e-2)
    end

    function test_step(s, d_set_values=zeros(3); dt=0.05, steps=5)
        s.integrator.ps[s.sys.stabilize] = true
        KiteModels.next_step!(s; dt=10.0)
        s.integrator.ps[s.sys.stabilize] = false
        @info "Stepping"
        for _ in 1:steps
            set_values = -s.set.drum_radius * s.integrator[s.sys.winch_force] + d_set_values
            KiteModels.next_step!(s, set_values; dt)
            @show s.integrator[s.sys.heading_x]
        end
    end

    @testset "Simulation Step with SysState" begin
        # Basic step and time advancement test
        KiteModels.reinit!(s, measure; prn=true, reload=false)
        sys_state_before = KiteModels.SysState(s)
        
        # Run a simulation step with zero set values
        set_values = [0.0, 0.0, 0.0]
        dt = 1/s.set.sample_freq
        t, _ = KiteModels.next_step!(s, set_values; dt=dt)
        KiteModels.update_sys_state!(sys_state_before, s)
        @test isapprox(t, dt, atol=TOL)
        
        # Run multiple steps
        num_steps = 10
        total_time = 0.0
        for _ in 1:num_steps
            step_time, _ = KiteModels.next_step!(s, set_values; dt=dt)
            total_time += step_time
        end
        sys_state_after = KiteModels.SysState(s)
        @test any(abs.(sys_state_after.X .- sys_state_before.X) .> 0.1)
        
        @testset "Course Direction at 70 Degrees Elevation" begin
            # Initialize at 70 degrees elevation
            measure.sphere_pos[1,:] .= deg2rad(70.0)
            @test measure.elevation ≈ deg2rad(70.0) atol=1e-6
            @test measure.azimuth ≈ 0.0 atol=1e-6
            
            KiteModels.reinit!(s, measure; prn=true)
            
            # Verify initial conditions
            @test s.integrator[s.sys.elevation] ≈ measure.elevation
            @test s.integrator[s.sys.azimuth] ≈ measure.azimuth
            
            # Run simulation steps
            test_step(s)
            
            # Check course direction
            sys_state = KiteModels.SysState(s)
            @info "Course at 70 deg elevation:" sys_state.course
            
            # At 70 degrees elevation, course should be roughly forward
            @show sys_state.course
            @test sys_state.course ≈ 0.0 atol=deg2rad(45.0)
        end

        # Utility function to calculate the signed angle difference between two angles
        function angle_diff(angle1, angle2)
            diff = angle1 - angle2
            # Normalize to [-π, π]
            while diff > π; diff -= 2π; end
            while diff < -π; diff += 2π; end
            return diff
        end
        
        @testset "Steering Response Using SysState" begin
            # Initialize model at moderate elevation
            measure.sphere_pos .= deg2rad.([70.0 70.0; 1.0 -1.0])
            KiteModels.reinit!(s, measure; prn=true, reload=false)
            test_step(s)
            sys_state_initial = KiteModels.SysState(s)
            
            # steering right
            KiteModels.reinit!(s, measure; prn=true, reload=false)
            test_step(s, [0, 5, -5]; steps=20)
            sys_state_right = KiteModels.SysState(s)
            
            # steering left
            KiteModels.reinit!(s, measure; prn=true, reload=false)
            test_step(s, [0, -5, 5]; steps=20)
            sys_state_left = KiteModels.SysState(s)
            
            # Check steering values
            @info "Steering:" sys_state_right.steering sys_state_left.steering
            @test sys_state_right.steering > 5.0
            @test sys_state_left.steering < -5.0
            
            # Check heading changes
            right_heading_diff = angle_diff(sys_state_right.heading, sys_state_initial.heading)
            @test right_heading_diff ≈ 2.0 atol=0.5
            left_heading_diff = angle_diff(sys_state_left.heading, sys_state_initial.heading)
            @test left_heading_diff ≈ -2.0 atol=0.5
        end
    end
end

# Restore original data path
set_data_path(old_path)

nothing
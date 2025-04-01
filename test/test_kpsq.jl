using Test, LinearAlgebra, KiteUtils, VortexStepMethod
using KiteModels
using Statistics

# Save original data path and set to package data directory
old_path = get_data_path()
set_data_path(joinpath(dirname(dirname(pathof(KiteModels))), "data"))

# Testing tolerance
const TOL = 1e-5
const BUILD_SYS = false

@testset verbose = true "RamAirKite MTK Model Tests" begin
    set = se("system_3l.yaml")
    wing = RamAirWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat"; 
                mass=set.mass, crease_frac=0.82, align_to_principal=true)
    aero = BodyAerodynamics([wing])
    vsm_solver = Solver(aero; solver_type=NONLIN, atol=1e-8, rtol=1e-8)
    

    # Utility functions for setup
    function create_test_model()
        set.segments = 2
        set.abs_tol = TOL
        set.rel_tol = TOL

        VortexStepMethod.init!(aero; init_aero=false)
        
        return RamAirKite(set, wing, aero, vsm_solver)
    end
    
    # @testset "Model Initialization Chain" begin
    #     s = create_test_model()
    #     if BUILD_SYS
    #         # Delete existing problem file to force init!
    #         prob_path = joinpath(get_data_path(), "prob.bin")
    #         if isfile(prob_path)
    #             @info "Removing existing serialized problem to test full initialization"
    #             rm(prob_path)
    #         end
            
    #         # 1. First time initialization - should create new model
    #         @info "Testing initial init! (should create new model)..."
    #         @time KiteModels.init!(s; prn=true)
            
    #         # Check that serialization worked
    #         @test isfile(prob_path)
            
    #         # Check initialization results
    #         @test !isnothing(s.integrator)
    #         @test !isnothing(s.sys)
    #         @test !isnothing(s.point_system)
    #     end
        
    #     # Keep references to first integrator and point system
    #     @test isnothing(s.integrator)
    #     first_integrator_ptr = objectid(s.integrator)
    #     first_point_system_ptr = objectid(s.point_system)
            
    #     # 2. First reinit! - should load from serialized file
    #     @info "Testing first reinit! (should load serialized file)..."
    #     @time KiteModels.reinit!(s; prn=true)
        
    #     # Check that it's a new integrator
    #     second_integrator_ptr = objectid(s.integrator)
    #     second_point_system_ptr = objectid(s.point_system)
    #     @test first_integrator_ptr != second_integrator_ptr
    #     @test first_point_system_ptr != second_point_system_ptr
            
    #     # 3. Second reinit! - should reuse existing integrator
    #     @info "Testing second reinit! (should reuse integrator)..."
    #     @time KiteModels.reinit!(s; prn=true)
        
    #     # This should create a new point system but reuse the existing integrator
    #     third_integrator_ptr = objectid(s.integrator)
    #     third_point_system_ptr = objectid(s.point_system)
    #     @test second_integrator_ptr == third_integrator_ptr # Should be the same 
    #     @test second_point_system_ptr == third_point_system_ptr # Should reuse the same object
            
    #     # Get positions from various sources
    #     pos_integrator, _, _, _, _, _, _, _, _, _, _ = s.get_state(s.integrator)
    #     sys_state = KiteModels.SysState(s)
        
    #     # Check dimension consistency
    #     @test size(pos_integrator, 2) == length(s.point_system.points)
    #     @test length(sys_state.X) == length(s.point_system.points)

        
    #     # Compare positions in different representations
    #     for (i, point) in enumerate(s.point_system.points)
    #         # Points' world positions should match integrator positions
    #         point_pos = point.pos_w
    #         integ_pos = pos_integrator[:, i]
    #         sys_state_pos = [sys_state.X[i], sys_state.Y[i], sys_state.Z[i]]
            
    #         @test isapprox(norm(point_pos), norm(integ_pos), rtol=1e-2)
    #         @test isapprox(norm(sys_state_pos), norm(integ_pos), rtol=1e-2)
            
    #         # Positions should not be zero (except ground points)
    #         if point.type != KiteModels.WINCH  # Skip ground points which might be at origin
    #             @test norm(point_pos) > 0.1
    #             @test norm(integ_pos) > 0.1
    #             @test norm(sys_state_pos) > 0.1
    #         end
    #     end
    # end
    
    # @testset "State Consistency" begin
    #     s = create_test_model()
    #     KiteModels.reinit!(s, prn=true)
    
    #     # Check quaternion normalization
    #     _, _, Q_b_w, elevation, azimuth, _, _, _, _, _, _ = s.get_state(s.integrator)
    #     @test isapprox(norm(Q_b_w), 1.0, atol=TOL)
    
    #     # Check elevation matches measurement
    #     @test isapprox(elevation, s.measure.elevation, atol=1e-2)
    
    #     # Change measurement and reinitialize
    #     old_elevation = s.measure.elevation
    #     s.measure.sphere_pos[1,:] .= deg2rad(85.0)
    #     KiteModels.reinit!(s, prn=true)
    
    #     # Get new state
    #     _, _, _, elevation_new, _, _, _, _, _, _, _ = s.get_state(s.integrator)
    
    #     # Verify state changed according to measurement
    #     @test !isapprox(elevation_new, old_elevation, atol=1e-2)
    #     @test isapprox(elevation_new, deg2rad(85.0), atol=1e-2)
    # end

    @testset "Simulation Step with SysState" begin
        # Basic step and time advancement test
        s = create_test_model()
        KiteModels.reinit!(s, prn=false)
        sys_state_before = KiteModels.SysState(s)
        
        # Run a simulation step with zero set values
        set_values = [0.0, 0.0, 0.0]
        dt = 1/s.set.sample_freq
        t, _ = KiteModels.next_step!(s; set_values=set_values, dt=dt)
        KiteModels.update_sys_state!(sys_state_before, s)
        @test isapprox(t, dt, atol=TOL)
        
        # Run multiple steps
        num_steps = 10
        total_time = 0.0
        for _ in 1:num_steps
            step_time, _ = KiteModels.next_step!(s; set_values=set_values, dt=dt)
            total_time += step_time
        end
        sys_state_after = KiteModels.SysState(s)
        @test any(abs.(sys_state_after.X .- sys_state_before.X) .> 0.1)
        
        @testset "Course Direction at Different Elevations" begin
            # Test at 80 degrees elevation
            s60 = create_test_model()
            s60.measure.sphere_pos[1,:] .= deg2rad(80.0)
            @test s60.measure.elevation ≈ deg2rad(80.0) atol=1e-6
            @test s60.measure.azimuth ≈ 0.0 atol=1e-6
            KiteModels.reinit!(s60; prn=false)
            @test s60.integrator[s60.sys.elevation] ≈ s60.measure.elevation
            @test s60.integrator[s60.sys.azimuth] ≈ s60.measure.azimuth
            for _ in 1:20
                set_values = -s60.set.drum_radius * s60.integrator[s60.sys.winch_force]
                KiteModels.next_step!(s60; set_values, dt)
                @show s60.integrator[s60.sys.kite_vel]
            end
            sys_state_60 = KiteModels.SysState(s60)
            @show sys_state_60.course s60.integrator[s60.sys.elevation_vel] s60.integrator[s60.sys.azimuth_vel] s60.integrator[s60.sys.kite_vel]
            @test sys_state_60.course ≈ 0.0 atol=0.2
            
            # Test at 90 degrees elevation
            s90 = create_test_model()
            s90.measure.sphere_pos[1,:] .= deg2rad(90.0)
            KiteModels.reinit!(s90, prn=false)
            for _ in 1:10
                set_values = -s90.set.drum_radius * s90.integrator[s90.sys.winch_force]
                KiteModels.next_step!(s90; set_values, dt)
            end
            sys_state_90 = KiteModels.SysState(s90)
            @info "Course at 90 deg:" sys_state_90.course
            @test abs(sys_state_90.course) ≈ π atol=0.2
        end
        
        # @testset "Steering Response Using SysState" begin
        #     # Initialize model at moderate elevation
        #     s_steer = create_test_model()
        #     s_steer.measure.sphere_pos[1,:] .= deg2rad(70.0)
        #     KiteModels.reinit!(s_steer, prn=false)
        #     for _ in 1:5
        #         KiteModels.next_step!(s_steer; set_values=[0.0, 0.0, 0.0], dt=dt)
        #     end
        #     sys_state_initial = KiteModels.SysState(s_steer)

        #     set_values_right = [0.0, 3.0, -3.0]  # Right steering
            
        #     # Create a copy of the model for steering right
        #     s_right = create_test_model()
        #     s_right.measure.sphere_pos[1,:] .= deg2rad(70.0)
        #     KiteModels.reinit!(s_right, prn=false)
            
        #     for _ in 1:5
        #         KiteModels.next_step!(s_right; set_values=[0.0, 0.0, 0.0], dt=dt)
        #     end
            
        #     for _ in 1:5
        #         KiteModels.next_step!(s_right; set_values=set_values_right, dt=dt)
        #     end
            
        #     # Get system state after steering right
        #     sys_state_right = KiteModels.SysState(s_right)
            
        #     # Create a copy of the model for steering left
        #     s_left = create_test_model()
        #     s_left.measure.sphere_pos[1,:] .= deg2rad(70.0)
        #     KiteModels.reinit!(s_left, prn=false)
            
        #     # Stabilize
        #     for _ in 1:5
        #         KiteModels.next_step!(s_left; set_values=[0.0, 0.0, 0.0], dt=dt)
        #     end
            
        #     # Test steering left (negative differential)
        #     set_values_left = [0.0, -3.0, 3.0]  # Left steering
            
        #     # Apply steering for several steps
        #     for _ in 1:5
        #         KiteModels.next_step!(s_left; set_values=set_values_left, dt=dt)
        #     end
            
        #     # Get system state after steering left
        #     sys_state_left = KiteModels.SysState(s_left)
            
        #     # Check steering values
        #     # Steering right should produce positive steering value
        #     @test sys_state_right.steering > 1.0
            
        #     # Steering left should produce negative steering value
        #     @test sys_state_left.steering < -1.0
            
        #     # Magnitudes should be approximately equal
        #     @test isapprox(abs(sys_state_right.steering), abs(sys_state_left.steering), rtol=0.3)
            
        #     # Check heading changes
        #     # After steering right, heading should change clockwise relative to initial
        #     right_heading_diff = angle_diff(sys_state_right.heading, sys_state_initial.heading)
        #     @test right_heading_diff > 0.05  # Should have some measurable change
            
        #     # After steering left, heading should change counterclockwise
        #     left_heading_diff = angle_diff(sys_state_left.heading, sys_state_initial.heading)
        #     @test left_heading_diff < -0.05  # Should have some measurable change
        # end
    end
    
    # Utility function to calculate the signed angle difference between two angles
    function angle_diff(angle1, angle2)
        diff = angle1 - angle2
        # Normalize to [-π, π]
        while diff > π; diff -= 2π; end
        while diff < -π; diff += 2π; end
        return diff
    end
    
    # @testset "System State Object" begin
    #     s = create_test_model()
    #     KiteModels.reinit!(s, prn=true)
    
    #     # Create system state object
    #     sys_state = KiteModels.SysState(s)
    
    #     # Check essential properties
    #     @test sys_state isa SysState
    #     @test length(sys_state.X) == length(s.point_system.points)
    #     @test length(sys_state.Y) == length(s.point_system.points)
    #     @test length(sys_state.Z) == length(s.point_system.points)
    #     @test length(sys_state.orient) == 4
    #     @test isapprox(norm(sys_state.orient), 1.0, atol=TOL)
    #     @test isfinite(sys_state.elevation)
    #     @test isfinite(sys_state.azimuth)
    #     @test isfinite(sys_state.heading)
    #     @test isfinite(sys_state.course)
    #     @test isfinite(sys_state.depower)
    #     @test isfinite(sys_state.steering)
    
    #     # Run a step and update system state
    #     KiteModels.next_step!(s; set_values=[0.0, 0.0, 0.0])
    #     KiteModels.update_sys_state!(sys_state, s)
    
    #     # Test updated values are reasonable
    #     @test all(isfinite, sys_state.X)
    #     @test all(isfinite, sys_state.Y)
    #     @test all(isfinite, sys_state.Z)
    #     @test all(isfinite, sys_state.vel_kite)
    #     @test isfinite(sys_state.v_reelout)
    
    #     # Compare with integrator positions
    #     pos_integrator, _, _, _, _, _, _, _, _, _, _ = s.get_state(s.integrator)
    
    #     for i in 1:length(s.point_system.points)
    #         integ_pos = pos_integrator[:, i]
    #         sys_state_pos = [sys_state.X[i], sys_state.Y[i], sys_state.Z[i]]
    
    #         # Check consistency between SysState and integrator
    #         @test isapprox(integ_pos[1], sys_state_pos[1], rtol=0.1)
    #         @test isapprox(integ_pos[2], sys_state_pos[2], rtol=0.1)
    #         @test isapprox(integ_pos[3], sys_state_pos[3], rtol=0.1)
    #     end
    # end
    
    # @testset "Aerodynamic Forces" begin
    #     s = create_test_model()
    #     KiteModels.reinit!(s, prn=true)
    
    #     # Run with VSM update every step
    #     set_values = [0.0, 0.0, 0.0]
    #     dt = 1/s.set.sample_freq
    #     KiteModels.next_step!(s; set_values=set_values, dt=dt, vsm_interval=1)
    
    #     # Get aerodynamic force and moment
    #     sys = s.sys
    #     aero_force = [s.integrator[sys.aero_force_b[i]] for i in 1:3]
    #     aero_moment = [s.integrator[sys.aero_moment_b[i]] for i in 1:3]
        
    #     # Check values are reasonable
    #     @test all(isfinite, aero_force)
    #     @test all(isfinite, aero_moment)
    #     @test norm(aero_force) > 0.1  # Should have some non-zero force
        
    #     # Test group moments
    #     group_moments = [s.integrator[sys.group_aero_moment[i]] for i in eachindex(s.point_system.groups)]
    #     @test all(isfinite, group_moments)
        
    #     # Test with different vsm_interval values
    #     KiteModels.reinit!(s, prn=true)
    #     KiteModels.next_step!(s; set_values=set_values, dt=dt, vsm_interval=5)
        
    #     # Force and moment should still be reasonable
    #     aero_force_2 = [s.integrator[sys.aero_force_b[i]] for i in 1:3]
    #     aero_moment_2 = [s.integrator[sys.aero_moment_b[i]] for i in 1:3]
        
    #     @test all(isfinite, aero_force_2)
    #     @test all(isfinite, aero_moment_2)
    # end
    
    # @testset "Control Response" begin
    #     s = create_test_model()
    #     KiteModels.reinit!(s, prn=true)
        
    #     # 1. Test power line control
    #     _, _, _, _, _, _, _, _, initial_tether_vel, _, _ = s.get_state(s.integrator)
        
    #     # Set non-zero reel-out speed
    #     set_values = [-10.0, 0.0, 0.0]  # First value controls main tether
    #     dt = 1/s.set.sample_freq
    #     num_steps = 10
        
    #     # Run with power line control
    #     for _ in 1:num_steps
    #         KiteModels.next_step!(s; set_values=set_values, dt=dt)
    #     end
        
    #     # Check tether velocity response
    #     _, _, _, _, _, _, _, _, final_tether_vel, _, _ = s.get_state(s.integrator)
    #     @test abs(final_tether_vel[1]) > abs(initial_tether_vel[1])
        
    #     # 2. Test steering control
    #     s = create_test_model()  # Fresh model
    #     KiteModels.reinit!(s, prn=true)
        
    #     # Get initial state
    #     _, _, _, _, _, _, initial_heading, _, _, initial_twist, _ = s.get_state(s.integrator)
        
    #     # Apply differential steering
    #     set_values = [0.0, 5.0, -5.0]  # Differential steering
        
    #     # Run with steering input
    #     for _ in 1:num_steps
    #         KiteModels.next_step!(s; set_values=set_values, dt=dt)
    #     end
        
    #     # Get heading and twist after steering
    #     _, _, _, _, _, _, final_heading, _, _, final_twist, _ = s.get_state(s.integrator)
        
    #     # Check twist angles differ between left and right
    #     @test !isapprox(final_twist[1], final_twist[end], atol=0.01)
        
    #     # Check twist differential increased from initial
    #     initial_twist_diff = abs(initial_twist[1] - initial_twist[end])
    #     final_twist_diff = abs(final_twist[1] - final_twist[end])
    #     @test final_twist_diff > initial_twist_diff
    # end
end

# Restore original data path
set_data_path(old_path)

nothing
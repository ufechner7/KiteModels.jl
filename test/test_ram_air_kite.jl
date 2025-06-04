# SPDX-FileCopyrightText: 2025 Uwe Fechner
#
# SPDX-License-Identifier: MIT

using Test, LinearAlgebra, KiteUtils, VortexStepMethod
using KiteModels
using Statistics

old_path = get_data_path()
package_data_path = joinpath(dirname(dirname(pathof(KiteModels))), "data")
temp_data_path = joinpath(tempdir(), "data")
Base.Filesystem.cptree(package_data_path, temp_data_path; force=true)
for bin_file in ["ram_air_kite_body_info.bin", "ram_air_kite_foil_polar.bin"]
    bin_path = joinpath(temp_data_path, bin_file)
    isfile(bin_path) && rm(bin_path)
end
set_data_path(temp_data_path)

# Testing tolerance
const TOL = 1e-5
const BUILD_SYS = true

@testset verbose = true "RamAirKite MTK Model Tests" begin
    set = se("system_ram.yaml")
    wing = RamAirWing(joinpath(get_data_path(), "ram_air_kite_body.obj"), joinpath(get_data_path(), "ram_air_kite_foil.dat"); 
                mass=set.mass, crease_frac=0.82, align_to_principal=true)
    aero = BodyAerodynamics([wing])
    vsm_solver = Solver(aero; solver_type=NONLIN, atol=1e-8, rtol=1e-8)
    point_system = PointMassSystem(set, wing)
    measure = Measurement()

    # Utility functions for setup
    function create_test_model()
        set.segments = 2
        set.abs_tol = TOL
        set.rel_tol = TOL
        VortexStepMethod.init!(aero; init_aero=false)
        return RamAirKite(set, aero, vsm_solver, point_system)
    end
    
    @testset "Model Initialization Chain" begin
        s = create_test_model()
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
        s = create_test_model()
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
        for _ in 1:steps
            set_values = -s.set.drum_radius * s.integrator[s.sys.winch_force] + d_set_values
            KiteModels.next_step!(s, set_values; dt)
        end
    end

    s = create_test_model()
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
        
        @testset "Course Direction at Different Elevations" begin
            # Test at 80 degrees elevation
            measure.sphere_pos[1,:] .= deg2rad(50.0)
            @test measure.elevation ≈ deg2rad(50.0) atol=1e-6
            @test measure.azimuth ≈ 0.0 atol=1e-6
            KiteModels.reinit!(s, measure; prn=true)
            @test s.integrator[s.sys.elevation] ≈ measure.elevation
            @test s.integrator[s.sys.azimuth] ≈ measure.azimuth
            test_step(s)
            sys_state_50 = KiteModels.SysState(s)
            @info "Course at 50 deg init elevation:" sys_state_50.course
            @test sys_state_50.course ≈ 0.0 atol=π/4
            
            # Test at 90 degrees elevation
            measure.sphere_pos[1,:] .= deg2rad(89.0)
            KiteModels.reinit!(s, measure; prn=true, reload=false)
            test_step(s)
            sys_state_89 = KiteModels.SysState(s)
            @info "Course at 89 deg init elevation:" sys_state_89.course
            @test abs(sys_state_89.course) ≈ π atol=π/4
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
            measure.sphere_pos[1,:] .= deg2rad(60.0)
            KiteModels.reinit!(s, measure; prn=true, reload=false)
            test_step(s)
            sys_state_initial = KiteModels.SysState(s)
            
            # steering right
            KiteModels.reinit!(s, measure; prn=true, reload=false)
            test_step(s, [0, 0, -5]; steps=20)
            sys_state_right = KiteModels.SysState(s)
            
            # steering left
            KiteModels.reinit!(s, measure; prn=true, reload=false)
            test_step(s, [0, -5, 0]; steps=20)
            sys_state_left = KiteModels.SysState(s)
            
            # Check steering values
            @info "Steering:" sys_state_right.steering sys_state_left.steering
            @test sys_state_right.steering > 2.0
            @test sys_state_left.steering < -2.0
            @test abs(sys_state_right.steering) ≈ abs(sys_state_left.steering) atol=1.0
            
            # Check heading changes
            right_heading_diff = angle_diff(sys_state_right.heading, sys_state_initial.heading)
            @test right_heading_diff > 0.6
            left_heading_diff = angle_diff(sys_state_left.heading, sys_state_initial.heading)
            @test left_heading_diff < -0.5
            @test abs(right_heading_diff) ≈ abs(left_heading_diff) atol=0.3
        end
    end
end

# Restore original data path
set_data_path(old_path)

nothing
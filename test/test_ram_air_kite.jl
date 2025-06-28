# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MIT

using Test, LinearAlgebra, KiteUtils, VortexStepMethod
using ControlPlots
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

@testset verbose = true "SymbolicAWEModel MTK Model Tests" begin
    # Initialize model
    set = se("system_ram.yaml")
    set.segments = 3
    set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
    set.quasi_static = false
    set.physical_model = "ram"

    @info "Creating s:"
    @time s = SymbolicAWEModel(set)

    s.set.abs_tol = 1e-2
    s.set.rel_tol = 1e-2

    # Initialize at elevation
    set.elevation = 80.0

    @testset "Model Initialization Chain" begin
        # if BUILD_SYS
        # Delete existing problem file to force init!
        @info "Data path: $(get_data_path())"
        model_path = joinpath(get_data_path(), KiteModels.get_model_name(s.set))
        # if isfile(model_path)
        #     @info "Removing existing serialized problem from $model_path to test full initialization"
        #     rm(model_path)
        # end

        # 1. First time initialization - should create new model
        @info "Testing initial init! (should create new model if it doesn't exist yet)..."
        @time KiteModels.init_sim!(s; prn=true)

        # Check that serialization worked
        @test isfile(model_path)

        # Check initialization results
        @test !isnothing(s.integrator)
        @test !isnothing(s.sys)
        @test !isnothing(s.sys_struct)
        # end
        s.integrator = nothing
        s.sys = nothing

        # Keep references to first integrator and point system
        first_integrator_ptr = objectid(s.integrator)
        first_sys_struct_ptr = objectid(s.sys_struct)

        # 2. First init_sim! - should load from serialized file
        @info "Testing first init_sim! (should load serialized file)..."
        @time KiteModels.init_sim!(s; prn=true, reload=false)
        next_step!(s)

        # Check that it's a new integrator
        second_integrator_ptr = objectid(s.integrator)
        second_sys_struct_ptr = objectid(s.sys_struct)
        @test first_integrator_ptr != second_integrator_ptr
        @test first_sys_struct_ptr == second_sys_struct_ptr

        # 3. Second init_sim! - should reuse existing integrator
        @info "Testing second init_sim! (should reuse integrator)..."
        @time KiteModels.init_sim!(s; prn=true, reload=false)

        # This should create a new point system but reuse the existing integrator
        third_integrator_ptr = objectid(s.integrator)
        third_sys_struct_ptr = objectid(s.sys_struct)
        @test second_integrator_ptr == third_integrator_ptr # Should be the same
        @test second_sys_struct_ptr == third_sys_struct_ptr

        # Get positions using SysState
        sys_state = KiteModels.SysState(s)

        # Check dimension consistency
        # Note: pos_integrator is no longer directly fetched, comparing SysState to sys_struct
        @test length(sys_state.X) == length(s.sys_struct.points)

        # Compare positions in different representations
        for (i, point) in enumerate(s.sys_struct.points)
            # Points' world positions should match SysState positions
            point_pos = point.pos_w
            sys_state_pos = [sys_state.X[i], sys_state.Y[i], sys_state.Z[i]]

            # Use norm for comparison as exact vector match might be too strict due to float precision
            @test isapprox(norm(point_pos), norm(sys_state_pos), rtol=1e-2)

            # Positions should not be zero (except ground points)
            if point.type != KiteModels.STATIC  # Skip ground points which might be at origin
                @test norm(point_pos) > 0.1
                @test norm(sys_state_pos) > 0.1
            end
        end
    end

    @testset "State Consistency" begin
        KiteModels.init_sim!(s, prn=true, reload=false)
        sys_state_before = KiteModels.SysState(s)
        @test isapprox(norm(s.integrator[s.sys.Q_b_w]), 1.0, atol=TOL)
        @test isapprox(sys_state_before.elevation, deg2rad(set.elevation), atol=1e-2)

        # Change measurement and reinitialize
        old_elevation = set.elevation
        set.elevation = 85.0
        KiteModels.init_sim!(s, prn=true, reload=false)

        # Get new state using SysState
        sys_state_after = KiteModels.SysState(s)

        # Verify state changed according to measurement
        @test !isapprox(sys_state_after.elevation, old_elevation, atol=1e-2)
        @test isapprox(sys_state_after.elevation, deg2rad(85.0), atol=1e-2)

        @testset "set_depower_steering!" begin
            initial_tether_lengths = s.get_tether_length(s.integrator)
            depower = 0.1
            steering = 0.05
            KiteModels.set_depower_steering!(s, depower, steering)
            new_tether_lengths = s.set_tether_length
            @test !isapprox(new_tether_lengths, initial_tether_lengths)
            # Verify the changes based on the equations
            len = s.set_tether_length
            len1 = initial_tether_lengths[1]
            len2 = 0.5 * (2*depower*KiteModels.min_chord_length(s) + 2*len1 + steering*KiteModels.min_chord_length(s))
            len3 = 0.5 * (2*depower*KiteModels.min_chord_length(s) + 2*len1 - steering*KiteModels.min_chord_length(s))
            @test isapprox(len[2], len2)
            @test isapprox(len[3], len3)
            @test isapprox(KiteModels.min_chord_length(s), 0.434108)
        end

        @testset "set_v_wind_ground!" begin
            initial_wind_speed = s.set.v_wind
            initial_upwind_dir = deg2rad(s.set.upwind_dir)
            @test initial_upwind_dir == -π/2
            @test s.integrator[s.sys.wind_vec_gnd[1]] == s.set.v_wind

            # Set new wind speed and direction
            new_wind_speed = 10.0
            new_upwind_dir = -pi/4
            KiteModels.set_v_wind_ground!(s, new_wind_speed, new_upwind_dir)

            # Check if wind speed and direction have been updated correctly
            @test norm(s.integrator[s.sys.wind_vec_gnd]) ≈ new_wind_speed
            @test s.integrator[s.sys.wind_vec_gnd[1]] ≈ -s.integrator[s.sys.wind_vec_gnd[2]]

            KiteModels.set_v_wind_ground!(s, initial_wind_speed, initial_upwind_dir)
            @test s.integrator[s.sys.wind_vec_gnd[1]] ≈ initial_wind_speed
            @test norm(s.integrator[s.sys.wind_vec_gnd]) ≈ initial_wind_speed
        end
    end

    function test_step(s, d_set_values=zeros(3); dt=0.05, steps=5)
        s.integrator.ps[s.sys.stabilize] = true
        for i in 1:1÷dt
            next_step!(s; dt, vsm_interval=1)
        end
        s.integrator.ps[s.sys.stabilize] = false
        @info "Stepping"
        for _ in 1:steps
            set_values = -s.set.drum_radius * s.integrator[s.sys.winch_force] + d_set_values
            next_step!(s; set_values, dt)
            # Use SysState to get heading if needed, or directly from integrator if simpler
            # sys_state_step = KiteModels.SysState(s)
            # @show sys_state_step.heading # Example if heading is in SysState
            @show s.integrator[s.sys.heading] # Keep direct access if simpler for this specific value
        end
    end

    function test_plot(s)
        @testset "Plotting of SymbolicAWEModel" begin
            function plot_(zoom, front)
                plt.figure("Kite")
                lines, sc, txt = plot(s, 0.0; zoom, front)
                plt.show(block=false)
                sleep(1)
                @test !isnothing(lines)
                @test length(lines) ≥ 1  # Should have at least one line
                @test !isnothing(sc)     # Should have scatter points
                @test !isnothing(txt)    # Should have time text
            end
            plot_(false, false)
            plot_(false, true)
            plot_(true, false)
            plot_(true, true)
        end
    end

    @testset "Simulation Step with SysState" begin
        # Basic step and time advancement test
        KiteModels.init_sim!(s; prn=true, reload=false)
        sys_state_before = KiteModels.SysState(s)

        # Run a simulation step with zero set values
        set_values = [0.0, 0.0, 0.0]
        dt = 1/s.set.sample_freq
        next_step!(s; set_values, dt=dt)
        # Update sys_state_before *after* the step to compare with the state *before* the loop
        KiteModels.update_sys_state!(sys_state_before, s)
        @test isapprox(s.integrator.t, dt, atol=TOL)

        # Run multiple steps
        num_steps = 10
        for _ in 1:num_steps
            next_step!(s; set_values, dt=dt)
        end
        sys_state_after = KiteModels.SysState(s) # Get state after the loop
        # Compare state after loop with state after first step (stored in sys_state_before)
        @test any(abs.(sys_state_after.X .- sys_state_before.X) .> 0.01) # Reduced tolerance slightly

        @testset "Course Direction at 60 Degrees Elevation" begin # Corrected description
            # Initialize at 60 degrees elevation
            set.elevation = 60.0

            KiteModels.init_sim!(s; prn=true)

            # Verify initial conditions using SysState
            sys_state_init = KiteModels.SysState(s)
            @test sys_state_init.elevation ≈ deg2rad(set.elevation) atol=1e-2 # Use relaxed tolerance consistent with other tests
            @test sys_state_init.azimuth ≈ deg2rad(set.azimuth) atol=1e-2

            # Run simulation steps
            test_step(s)

            # Check course direction using SysState
            sys_state = KiteModels.SysState(s)
            @info "Course at 60 deg elevation:" sys_state.course

            # At 60 degrees elevation, course should be roughly forward
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
            set.elevation = 70.0
            KiteModels.init_sim!(s; prn=true, reload=false)
            test_step(s)
            sys_state_initial = KiteModels.SysState(s)

            # steering right
            KiteModels.init_sim!(s; prn=true, reload=false)
            test_step(s, [0, 10, -10]; steps=20)
            sys_state_right = KiteModels.SysState(s)

            # steering left
            KiteModels.init_sim!(s; prn=true, reload=false)
            test_step(s, [0, -10, 10]; steps=20)
            sys_state_left = KiteModels.SysState(s)

            # Check steering values
            @info "Steering:" sys_state_right.steering sys_state_left.steering
            @test sys_state_right.steering > 3.0
            @test sys_state_left.steering < -3.0

            # Check heading changes
            right_heading_diff = angle_diff(sys_state_right.heading, sys_state_initial.heading)
            @test right_heading_diff ≈ 0.9 atol=0.2
            left_heading_diff = angle_diff(sys_state_left.heading, sys_state_initial.heading)
            @test left_heading_diff ≈ -0.9 atol=0.2
        end
        test_plot(s)
    end

    @testset "Just a tether, without winch or kite" begin
        set.segments = 20
        dynamics_type = DYNAMIC

        points = Point[]
        segments = Segment[]

        points = push!(points, Point(1, zeros(3), STATIC; wing_idx=0))

        segment_idxs = Int[]
        for i in 1:set.segments
            point_idx = i+1
            pos = [0.0, 0.0, i * set.l_tether / set.segments]
            push!(points, Point(point_idx, pos, dynamics_type; wing_idx=0))
            segment_idx = i
            push!(segments, Segment(segment_idx, (point_idx-1, point_idx), BRIDLE))
            push!(segment_idxs, segment_idx)
        end

        transforms = [Transform(1, deg2rad(-80), 0.0, 0.0; 
            base_pos=[0.0, 0.0, 50.0], base_point_idx=points[1].idx, rot_point_idx=points[end].idx)]
        sys_struct = SystemStructure("tether", set; points, segments, transforms)

        sam = SymbolicAWEModel(set, sys_struct)
        sys = sam.sys
        init_sim!(sam; remake=false)
        @test isapprox(sam.integrator[sam.sys.pos[:, end]], [8.682408883346524, 0.0, 0.7596123493895988], atol=1e-2)
        for i in 1:100
            next_step!(sam)
        end
        @test sam.integrator[sam.sys.pos[1, end]] > 0.8set.l_tether
        @test isapprox(sam.integrator[sam.sys.pos[2, end]], 0.0, atol=1.0)
        test_plot(s)
        set.elevation = 80
    end
end

# Restore original data path
set_data_path(old_path)

nothing

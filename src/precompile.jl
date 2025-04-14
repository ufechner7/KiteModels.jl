function decompress_binary(infile, outfile; chunksize=4096)
    open(infile) do input
        open(outfile, "w") do output
            stream = XzDecompressorStream(input)
            while !eof(stream)
                write(output, read(stream, chunksize))
            end
        end
    end
end

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    path = dirname(pathof(@__MODULE__))
    set_data_path(joinpath(path, "..", "data"))

    set = se("system.yaml")
    set.kcu_diameter = 0
    kps4_::KPS4 = KPS4(KCU(set))
    kps3_::KPS3 = KPS3(KCU(se("system.yaml")))

    input_file = joinpath(path, "..", "data", "prob_dynamic_3_seg.bin.default.xz")
    output_file = joinpath(path, "..", "data", "prob_dynamic_3_seg.bin.default")
    if isfile(input_file) && ! isfile(output_file)
        using CodecXz
        decompress_binary(input_file, output_file)
    end

    @assert ! isnothing(kps4_.wm)
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether(
        # they belong to your package or not (on Julia 1.8 and higher)
        integrator = KiteModels.init_sim!(kps3_; stiffness_factor=0.035, prn=false)
        integrator = KiteModels.init_sim!(kps4_; delta=0.03, stiffness_factor=0.05, prn=false)

        if isfile(joinpath(path, "..", "data", "prob_dynamic_3_seg.bin.default"))
            local set
            # Simulation parameters
            dt = 0.05
            total_time = 10  # Longer simulation to see oscillations
            vsm_interval = 2
            steps = Int(round(total_time / dt))

            # Steering parameters
            steering_freq = 1/2  # Hz - full left-right cycle frequency
            steering_magnitude = 5.0      # Magnitude of steering input [Nm]

            # Initialize model
            set = se("system_ram.yaml")
            set.segments = 3
            set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
            set.quasi_static = false

            wing = RamAirWing(set; prn=false)
            aero = BodyAerodynamics([wing])
            vsm_solver = Solver(aero; solver_type=NONLIN, atol=2e-8, rtol=2e-8)
            point_system = generate_ram_point_system(set, wing)
            s = RamAirKite(set, aero, vsm_solver, point_system)

            measure = Measurement()
            s.set.abs_tol = 1e-5
            s.set.rel_tol = 1e-4

            # Initialize at elevation
            measure.sphere_pos .= deg2rad.([60.0 60.0; 1.0 -1.0])
            KiteModels.init_sim!(s, measure; prn=false, precompile=true)
        end
        nothing
    end
end
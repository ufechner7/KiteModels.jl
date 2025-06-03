# SPDX-License-Identifier: MIT

#= MIT License

Copyright (c) 2024 Uwe Fechner, Bart van de Lint

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. =#

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

function filecmp(path1::AbstractString, path2::AbstractString)
    stat1, stat2 = stat(path1), stat(path2)
    if !(isfile(stat1) && isfile(stat2)) || filesize(stat1) != filesize(stat2)
        return false
    end
    open(path1, "r") do file1
        open(path2, "r") do file2
            buf1 = Vector{UInt8}(undef, 32768)
            buf2 = similar(buf1)
            while !eof(file1) && !eof(file2)
                n1 = readbytes!(file1, buf1)
                n2 = readbytes!(file2, buf2)
                if n1 != n2 || buf1[1:n1] != buf2[1:n2]
                    return false
                end
            end
            return eof(file1) && eof(file2)
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
    ver = "$(VERSION.major).$(VERSION.minor)_"

    input_file  = joinpath(path, "..", "data", "prob_dynamic_" * ver * "3_seg.bin.default.xz")
    prob_file   = joinpath(path, "..", "data", "prob_dynamic_" * ver * "3_seg.bin")
    output_file = joinpath(prob_file * ".default")
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
        if VERSION.minor == 11
            m1 = "Manifest-v1.11.toml"
            m2 = "Manifest-v1.11.toml.default"
        else
            m1 = "Manifest-v1.10.toml"
            m2 = "Manifest-v1.10.toml.default"
        end
        if isfile(output_file) && filecmp(m1, m2)
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
            point_system = PointMassSystem(set, wing)
            s = RamAirKite(set, aero, vsm_solver, point_system)

            measure = Measurement()
            s.set.abs_tol = 1e-5
            s.set.rel_tol = 1e-4

            # Initialize at elevation
            measure.sphere_pos .= deg2rad.([60.0 60.0; 1.0 -1.0])
            KiteModels.init_sim!(s, measure; prn=false, precompile=true)
            cp(output_file, prob_file; force=true)
        end
        nothing
    end
end
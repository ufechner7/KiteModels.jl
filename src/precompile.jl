# Copyright (c) 2024 Uwe Fechner, Bart van de Lint
# SPDX-License-Identifier: MIT

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

    set = load_settings("system.yaml")
    set.kcu_diameter = 0.0
    kps4_::KPS4 = KPS4(KCU(set))
    kps3_::KPS3 = KPS3(KCU(set))
    ver = "$(VERSION.major).$(VERSION.minor)_"

    @assert ! isnothing(kps4_.wm)
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether(
        # they belong to your package or not (on Julia 1.8 and higher)
        integrator = KiteModels.init_sim!(kps3_; stiffness_factor=0.035, prn=false)
        integrator = KiteModels.init_sim!(kps4_; delta=0.03, stiffness_factor=0.05, prn=false)

        sam_set = load_settings("system_ram.yaml")
        sam_set.segments = 3
        set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
        sam_set.quasi_static = false
        sam_set.physical_model = "ram"
        model_name   = get_model_name(sam_set)
        model_file   = normpath(joinpath(path, "..", "data", model_name))
        output_file = normpath(joinpath(path, "..", "data", model_name * ".default"))
        input_file  = normpath(joinpath(path, "..", "data", model_name * ".default.xz"))
        if isfile(input_file) && ! isfile(output_file)
            using CodecXz
            decompress_binary(input_file, output_file)
            @info "Decompressed $input_file to $output_file"
        elseif isfile(output_file)
            @info "Output file $output_file already exists, skipping decompression."
        else
            @error "Input file $input_file does not exist, skipping decompression."
        end
        if VERSION.minor == 11
            m1 = "Manifest-v1.11.toml"
            m2 = "Manifest-v1.11.toml.default"
        else
            m1 = "Manifest-v1.10.toml"
            m2 = "Manifest-v1.10.toml.default"
        end
        if filecmp(m1, m2)
            @info "Manifest files match, using the default xz files will work!"
        else
            @warn "Manifest files differ, no precompilation will be done."
        end
        # Check if the output file exists and is the same as the input file
        if isfile(output_file) && filecmp(m1, m2)
            s = SymbolicAWEModel(sam_set)

            # Initialize at elevation
            KiteModels.init_sim!(s; prn=false, precompile=true)
            @info "Copying $output_file to $model_file !"
            cp(output_file, model_file; force=true)
            find_steady_state!(s)
            steps = Int(round(10 / 0.05))
            logger = Logger(length(s.sys_struct.points), steps)
            sys_state = SysState(s)
        end
        nothing
    end
end   
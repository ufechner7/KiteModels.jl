using KiteModels, KiteUtils
using Test

cd("..")
KiteUtils.set_data_path("") 
@testset verbose = true "Testing KiteModels..." begin
    include("test_orientation.jl")
    include("test_kps3.jl")
    include("test_kps4.jl")
    # include("test_ram_air_kite.jl")
    include("bench3.jl")
    include("bench4.jl")
    include("test_helpers.jl")
    include("test_inertia_calculation.jl")
end
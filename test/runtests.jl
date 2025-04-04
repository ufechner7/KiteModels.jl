using KiteModels, KiteUtils
using Test

cd("..")
KiteUtils.set_data_path("") 
@testset verbose = true "Testing KiteModels..." begin
    include("test_orientation.jl")
    include("test_kps3.jl")
    include("test_kps4.jl")
    if ! haskey(ENV, "NO_MTK")  
        include("test_kps4_3l.jl")
    end
    include("bench3.jl")
    include("bench4.jl")
    include("test_helpers.jl")
    include("test_inertia_calculation.jl")
end
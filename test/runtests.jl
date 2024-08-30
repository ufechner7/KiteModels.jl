using KiteModels, KiteUtils
using Test

cd("..")
KiteUtils.set_data_path("") 
include("test_kps3.jl")
include("test_kps4.jl")
include("test_kps4_3l.jl")
include("bench3.jl")
include("bench4.jl")
include("test_helpers.jl")

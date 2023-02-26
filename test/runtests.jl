using KiteModels, KiteUtils
using Test

cd("..")
KiteUtils.set_data_path("") 
include("test_kps3.jl")
# include("test_kps4.jl")
include("bench3.jl")
# include("bench4.jl")

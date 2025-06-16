# SPDX-FileCopyrightText: 2025 Uwe Fechner
#
# SPDX-License-Identifier: MIT

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using KiteModels, StaticArrays, LinearAlgebra, BenchmarkTools
function test(vec)
    KiteModels.norm(vec)
end
vec=MVector{3}([1.0,2,3])

test(vec)
bytes=@allocated test(vec)
println("Allocate $bytes bytes!")
@benchmark test($vec)
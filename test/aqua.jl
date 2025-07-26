# SPDX-FileCopyrightText: 2024 Uwe Fechner
# SPDX-License-Identifier: MIT

using Pkg
if ! ("Aqua" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end

using KiteModels, Aqua, Test
@testset "Aqua.jl" begin
    Aqua.test_all(
      KiteModels;
      stale_deps=(ignore=[:PyCall, :CodecXz, :REPL],), # CodecXz is used during precompilation only
      deps_compat=(ignore=[:PyCall],),                 # PyCall is needed for CI to recompile Python
      piracies=false                                   # the norm function is doing piracy for performance reasons
    )
end

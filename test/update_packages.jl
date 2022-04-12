@info "Updating packages ..."
using Pkg
Pkg.instantiate()
Pkg.update()
Pkg.precompile()

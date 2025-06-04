# SPDX-FileCopyrightText: 2025 Uwe Fechner
#
# SPDX-License-Identifier: MIT

@info "Updating packages ..."
using Pkg
Pkg.instantiate()
Pkg.update()
Pkg.precompile()

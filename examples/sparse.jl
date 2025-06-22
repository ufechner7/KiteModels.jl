# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

using DifferentiationInterface
using SparseArrays
using ADTypes

# 1. Define your function
function f(x)
    y = similar(x, length(x)-1)
    for i in 1:length(y)
        y[i] = x[i+1]^2 - x[i]^2
    end
    return y
end

# 2. Choose a finite difference backend
backend = AutoFiniteDiff()

# 3. Provide the known sparsity pattern as a SparseMatrixCSC{Bool}
#    For f: ℝⁿ → ℝⁿ⁻¹, each output depends only on x[i] and x[i+1]
n = 5
rows = Int[]
cols = Int[]
for i in 1:n-1
    push!(rows, i); push!(cols, i)     # y[i] depends on x[i]
    push!(rows, i); push!(cols, i+1)   # y[i] depends on x[i+1]
end
S = sparse(rows, cols, trues(length(rows)), n-1, n) # Bool-valued sparsity pattern

# 4. Prepare the Jacobian using the known sparsity pattern
jac_prep = prepare_jacobian(f, backend, zeros(n); sparsity_pattern=S)

# 5. Compute the sparse Jacobian at a point
x = randn(n)
J = jacobian(f, jac_prep, backend, x)

# 6. Show the result
println("Jacobian (sparse):")
display(J)

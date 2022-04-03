using NOMAD, LinearAlgebra

function f(x)
  return Float64[x[1]^2 + x[2]^2, (x[1]-2)^2]
end

function eval_fct(x)
  bb_outputs = [norm(f(x))]
  success = true
  count_eval = true
  return (success, count_eval, bb_outputs)
end

pb = NomadProblem(2, # number of inputs of the blackbox
                  1, # number of outputs of the blackbox
                  ["OBJ"], # type of outputs of the blackbox
                  eval_fct;
                  lower_bound=[-5.0, -5.0],
                  upper_bound=[5.0, 5.0])

pb.options.display_stats = ["BBE", "OBJ", "SOL"] 

result = solve(pb, [3.0, 3.0])
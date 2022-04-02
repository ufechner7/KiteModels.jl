using NOMAD

function f(x)
  return x[1]^2 + x[2]^2
end

function c(x)
  return 1 - x[1]
end

function eval_fct(x)
  bb_outputs = [f(x), c(x)]
  success = true
  count_eval = true
  return (success, count_eval, bb_outputs)
end

pb = NomadProblem(2, # number of inputs of the blackbox
                  2, # number of outputs of the blackbox
                  ["OBJ", "EB"], # type of outputs of the blackbox
                  eval_fct;
                  lower_bound=[-5.0, -5.0],
                  upper_bound=[5.0, 5.0])

result = solve(pb, [3.0, 3.0])
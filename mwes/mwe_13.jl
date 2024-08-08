using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

function decay(; name)
    @parameters a
    @variables x(t) f(t)
    ODESystem([
            D(x) ~ -a * x + f
        ], t;
        name = name)
end

@named decay1 = decay()
@named decay2 = decay()

connected = compose(
    ODESystem([decay2.f ~ decay1.x
               D(decay1.f) ~ 0], t; name = :connected), decay1, decay2)

equations(connected)
simplified_sys = structural_simplify(connected)
equations(simplified_sys)

x0 = [decay1.x => 1.0
      decay1.f => 0.0
      decay2.x => 1.0]
p = [decay1.a => 0.1
     decay2.a => 0.2]

using OrdinaryDiffEq
prob = ODEProblem(simplified_sys, x0, (0.0, 100.0), p)
sol = solve(prob, Tsit5())
@time sol[decay2.f][end]
@time sol.u[end][2]
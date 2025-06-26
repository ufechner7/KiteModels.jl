using ModelingToolkit, Dierckx, Parameters

function dirac_delta(t, t0)
    return ifelse(t == t0, t, 0)
end

function step(t, t0)
    return ifelse(t >= t0, t, 0)
end

@register_symbolic dirac_delta(t, t0)
@register_symbolic step(t, t0)

@with_kw mutable struct TestKPS
    "Function for calculation the lift coefficent, using a spline based on the provided value pairs."
    calc_cl::Spline1D
    "Function for calculation the drag coefficent, using a spline based on the provided value pairs."
    calc_cd::Spline1D
end
function KPS_TEST(cl_alpha, cl_list, cd_alpha, cd_list)
    s = TestKPS(calc_cl=Spline1D(cl_alpha, cl_list), calc_cd=Spline1D(cd_alpha, cd_list))
    return s
end

function Cl(s::TestKPS, alpha)
    return s.calc_cl(alpha)
end

function Cd(s::TestKPS, alpha)
    return s.calc_cd(alpha)
end

@register_symbolic Cl(s::TestKPS, alpha)
@register_symbolic Cd(s::TestKPS, alpha)

function model(s::TestKPS)
    @variables alpha1p(ModelingToolkit.t_nounits)[1:1]
    cl = Cl(s, alpha1p[1])
    cd = Cd(s, alpha1p[1])
    return cl, cd
end

KPS = KPS_TEST(
    [-5.0,   0.0,  5.0,  7.0],
    [ 0.0,   0.0,  0.08, 0.125],
    [-20.0,  0.0,  7.5, 8.0],
    [  0.5,  1.0,  0.06, 0.07]
)

println(model(KPS))
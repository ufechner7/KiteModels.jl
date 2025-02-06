#= MIT License

Copyright (c) 2024 Uwe Fechner and Bart van de Lint

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. =#

#= Model of a kite-power system in implicit form: residual = f(y, yd)

This model implements a 3D mass-spring system with reel-out. It uses six tether segments (the number can be
configured in the file data/settings.yaml). The kite is modelled using 4 point masses and 2n aerodynamic 
surfaces. The spring constant and the damping decrease with the segment length. The aerodynamic kite forces
are acting on the inertial kite point. 

Four point kite model, included from KiteModels.jl.

Scientific background: http://arxiv.org/abs/1406.6218 =#

# Array of connections of bridlepoints.
# First point, second point, unstressed length.
const SPRINGS_INPUT_3L = [1.      4.  -1. # s1: E, A
                        1.      2.  -1. # s2, E, C
                        1.      3.  -1. # s3, E, D
                        2.      3.  -1. # s4, C, D
                        2.      4.  -1. # s5, C, A
                        3.      4.  -1. # s6, D, A
                        ]
# E = 1, C = 2, D = 3, A = 4
# E = segments*3+1, C = segments*3+2, D = segments*3+3, A = segments*3+4

# TODO: add line multiplier: multiple lines doing same thing
const KITE_SPRINGS_3L = 6
const KITE_PARTICLES_3L = 4
const MeasureFloat = Float32

@with_kw mutable struct Measurements
    set_values::MVector{3, MeasureFloat}    = [-1., -1., -50.]
    tether_length::MVector{3, MeasureFloat} = [51., 51., 49.]
    tether_vel::MVector{3, MeasureFloat}    = zeros(3)
    tether_acc::MVector{3, MeasureFloat}    = zeros(3)
    "elevation and azimuth in spherical coordinate system with columns (left, right) and rows (elevation, azimuth)"
    sphere_pos::Matrix{SimFloat}            = deg2rad.([89.0 89.0; 1.0 -1.0])
    sphere_vel::Matrix{SimFloat}            = zeros(SimFloat, 2, 2)
    sphere_acc::Matrix{SimFloat}            = zeros(SimFloat, 2, 2)
end

"""
    mutable struct KPSQ{S, T, P, Q, SP} <: AbstractKiteModel

State of the kite power system, using a quaternion kite model and three steering lines to the ground. Parameters:
- S: Scalar type, e.g. SimFloat
  In the documentation mentioned as Any, but when used in this module it is always SimFloat and not Any.
- V: Vector type, e.g. KVec3
- P: number of tether points of the system, 3segments+3
Normally a user of this package will not have to access any of the members of this type directly,
use the input and output functions instead.

$(TYPEDFIELDS)
"""
@with_kw mutable struct KPSQ{S, V, P} <: AbstractKiteModel # TODO: subdivide in kite-changing fields and non-kite-changing fields. combine this with fast caching
    "Reference to the settings struct"
    set::Settings
    "The last initial elevation"
    last_init_elevation::S     = 0.0
    "The last initial tether length"
    last_init_tether_length::S = 0.0
    "Reference to the last settings hash"
    last_set_hash::UInt64   = 0
    "Reference to the last measurement hash"
    last_measure_hash::UInt64 = 0
    "Reference to the atmospheric model as implemented in the package AtmosphericModels"
    am::AtmosphericModel = AtmosphericModel()
    "Function for calculating the lift coefficent, using linear interpolation based on the provided value pairs."
    cl_interp::Function
    "Function for calculating the drag coefficent, using linear interpolation based on the provided value pairs."
    cd_interp::Function
    "Function for calculating the trailing edge force coefficient, using linear interpolation based on the provided value pairs."
    c_te_interp::Function
    "Reference to the motor models as implemented in the package WinchModels. index 1: middle motor, index 2: left motor, index 3: right motor"
    motors::SizedArray{Tuple{3}, AbstractWinchModel}
    "wind vector at the height of the kite" 
    v_wind::V =           zeros(S, 3)
    "wind vector at reference height" 
    v_wind_gnd::V =       zeros(S, 3)
    "wind vector used for the calculation of the tether drag"
    v_wind_tether::V =    zeros(S, 3)
    "apparent wind vector at the kite"
    v_apparent::V =       zeros(S, 3)
    "tether positions"
    pos::Matrix{S} = zeros(S, 3, P)
    "unstressed segment lengths of the three tethers [m]"
    segment_lengths::V =           zeros(S, 3)
    "relative start time of the current time interval"
    t_0::S =               0.0
    "unstretched tether length"
    tether_lengths::V =          zeros(S, 3)
    "air density at the height of the kite"
    rho::S =               0.0
    "tether masses"
    masses::V         = zeros(S, P)
    "unit spring coefficient"
    c_spring::V = zeros(S, 3)
    "unit damping coefficient"
    damping::V = zeros(S, 3)
    "whether or not to use torque control instead of speed control"
    torque_control::Bool = false
    "x vector of kite reference frame"
    e_x::V =                 zeros(S, 3)
    "y vector of kite reference frame"
    e_y::V =                 zeros(S, 3)
    "z vector of kite reference frame"
    e_z::V =                 zeros(S, 3)
    "Point index of A - trailing edge point"
    i_A::Int64 =           0
    "Point index of B - trailing edge point"
    i_B::Int64 =           0
    "Point index of C - middle tether last point"
    i_C::Int64 =           0
    "Angle of left tip"
    γ_l::S =     0.0
    "Angle of point C"
    γ_D::S =     0.0
    "Kite length at point C"
    kite_length_D::S =     0.0
    "Simplified system of the mtk model"
    simple_sys::Union{ModelingToolkit.ODESystem, Nothing} = nothing
    "Velocity of the kite"
    vel_kite::V =           zeros(S, 3)
    "Initial torque or speed set values"
    init_set_values::V =    zeros(S, 3)
    "Smooth sign constant"
    ϵ::S =      0.0
    "Relative damping for flaps"
    flap_damping::S     = 0.75
    "Measured data points used to create an initial state"
    measure::Measurements = Measurements()
    "Buffer for expected kite pos and vel, defined as (pos, vel), (x, y, z), (value)"
    expected_tether_pos_vel_buffer::Array{S, 3} = zeros(S, 2, 3, (P ÷ 3))
    "Buffers for Jacobian for expected tether point velocities relative to winch velocities and kite velocities"
    J_buffer::Matrix{S} = zeros(S, P, 2)
    "Makes autodiff faster"
    prep::Union{Nothing, Any} = nothing
    "Buffer for jacobian y values (vectorized velocities)"
    y_buffer::V = zeros(S, P)
    "Buffer for jacobian x values (angle, distance)"
    x_buffer::V = zeros(S, 2)
    "Inertia around kite x y and z axis of the principal frame"
    I_p::V = zeros(S, 3)
    "Inertia around kite x y and z axis of the body frame"
    I_b::V = zeros(S, 3)
    "Damping of the kite rotation"
    orient_damping::S = zero(S)
    "rotation from kite body frame to world frame"
    Q_p_w::V = zeros(S, 4)
    "rotation from kite principal frame to body frame"
    Q_p_b::V = zeros(S, 4)
    "rotation from kite body frame to kite principal frame"
    R_b_p::Matrix{S} = zeros(S, 3, 3)
    "translation from kite point to circle center in body frame along positive z"
    pos_circle_center_b::V = zeros(S, 3)
    "center of mass position of the kite"
    kite_pos::V = zeros(S, 3)
    "quaternion orientation of the kite"
    q::V = zeros(S, 4)
    "aero points center of mass positions in principal frame"
    seg_com_pos_p::Matrix{S} = zeros(S, 3, 2set.aero_surfaces)
    "aero points center of pressure positions in principal frame"
    seg_cop_pos_p::Matrix{S} = zeros(S, 3, 2set.aero_surfaces)
    "cop pos in body frame"
    seg_cop_pos_b::Matrix{S} = zeros(S, 3, 2set.aero_surfaces)
    "last left tether point position in body frame"
    pos_A_b::V = zeros(S, 3)
    "last right tether point position in body frame"
    pos_B_b::V = zeros(S, 3)
    "last middle tether point position in body frame"
    pos_C_b::V = zeros(S, 3)
    "last middle tether point position in principal frame"
    pos_C_p::V = zeros(S, 3)
    "mass of each kite segment"
    seg_mass::V = zeros(S, 2set.aero_surfaces)
    "A point projected onto kite in z-axis in body frame"
    pos_D_b::V = zeros(S, 3)
    "B point projected onto kite in z-axis in body frame"
    pos_E_b::V = zeros(S, 3)
    "Initialization values for kite state"
    u0map::Union{Vector{Pair{Num, S}}, Nothing} = nothing
    "Initialization values for kite parameters"
    p0map::Union{Vector{Pair{Num, S}}, Nothing} = nothing
    "Distance of the kite com from winch"
    distance::S = zero(S)
    "Angle of the trailing edges of the kite"
    te_angle::V = zeros(S, 2)

    set_Q_p_w::Function             = (val, prob) -> nothing
    set_ω_p::Function               = (val, prob) -> nothing
    set_kite_pos::Function          = (val, prob) -> nothing 
    set_kite_vel::Function          = (val, prob) -> nothing
    set_pos::Function               = (val, prob) -> nothing
    set_vel::Function               = (val, prob) -> nothing
    set_trailing_edge_angle::Function = (val, prob) -> nothing
    set_trailing_edge_ω::Function    = (val, prob) -> nothing
    set_gust_factor::Function        = (val, prob) -> nothing
    set_tether_length::Function      = (val, prob) -> nothing
    set_tether_vel::Function         = (val, prob) -> nothing
    set_v_wind_gnd::Function        = (val) -> nothing
    set_v_wind::Function            = (val) -> nothing
    set_set_values::Function        = (val) -> nothing

    get_distance_acc::Function       = () -> nothing
    get_pos::Function               = () -> nothing
    get_vel::Function               = () -> nothing
    get_acc::Function               = () -> nothing
    get_trailing_edge_angle::Function = () -> nothing
    get_trailing_edge_α::Function    = () -> nothing
    get_kite_pos::Function          = () -> nothing
    get_kite_vel::Function          = () -> nothing
    get_kite_acc::Function          = () -> nothing
    get_R_b_w::Function             = () -> nothing
    get_Q_p_w::Function             = () -> nothing
    get_e_x::Function               = () -> nothing
    get_e_y::Function               = () -> nothing
    get_e_z::Function               = () -> nothing
    get_tether_force::Function      = () -> nothing
    get_tether_length::Function     = () -> nothing
    get_tether_vel::Function        = () -> nothing
    get_tether_acc::Function        = () -> nothing
    get_kite_force::Function        = () -> nothing
    get_kite_torque_p::Function     = () -> nothing
    get_heading::Function           = () -> nothing
    get_ω_b::Function              = () -> nothing
    get_α_b::Function              = () -> nothing
    get_force::Function            = () -> nothing
    get_e_te_A::Function           = () -> nothing
    get_e_te_B::Function           = () -> nothing
    get_elevation_vel::Function     = () -> nothing
    get_elevation_acc::Function     = () -> nothing
    get_azimuth_vel::Function      = () -> nothing
    get_azimuth_acc::Function      = () -> nothing

    prob::Union{OrdinaryDiffEqCore.ODEProblem, Nothing} = nothing
    integrator::Union{OrdinaryDiffEqCore.ODEIntegrator, Sundials.CVODEIntegrator, Nothing} = nothing
end

function kite_torque_b(s)
    return rotate_by_quaternion(kite_torque_p(s), s.Q_p_b)
end

"""
    clear!(s::KPSQ)

Initialize the kite power model.
"""
function clear!(s::KPSQ)
    P = 3s.set.segments + 3
    S = eltype(s.pos)
    s.pos = zeros(S, 3, P)
    s.masses = zeros(S, P)
    s.expected_tether_pos_vel_buffer = zeros(S, 2, 3, (P ÷ 3))
    s.J_buffer = zeros(S, P, 2)
    s.y_buffer = zeros(S, P)
    
    s.t_0 = 0.0                              # relative start time of the current time interval
    # s.last_reel_out_speeds = zeros(3)
    s.v_wind_gnd    .= [s.set.v_wind, 0.0, 0.0]    # wind vector at reference height
    s.v_wind_tether .= [s.set.v_wind, 0.0, 0.0]
    s.v_apparent    .= [s.set.v_wind, 0.0, 0.0]
    height = sin(deg2rad(s.set.elevation)) * (s.set.l_tether)
    s.v_wind .= s.v_wind_gnd * calc_wind_factor(s.am, height)
    s.e_x .= 0.0
    s.e_y .= 0.0
    s.e_z .= 0.0
    s.tether_lengths .= [s.set.l_tether for _ in 1:3]
    s.γ_l = π/2 - s.set.min_steering_line_distance/(2*s.set.radius)
    s.segment_lengths .= s.tether_lengths ./ s.set.segments
    s.i_A = s.set.segments*3+1
    s.i_B = s.set.segments*3+2
    s.i_C = s.set.segments*3+3
    s.rho = s.set.rho_0
    c_spring = s.set.e_tether * (s.set.d_tether/2000.0)^2 * pi
    s.c_spring .= [c_spring, c_spring, 2c_spring]
    s.damping .= (s.set.damping / s.set.c_spring) * s.c_spring
    init_masses!(s)

    width, radius, tip_length, middle_length = s.set.width, s.set.radius, s.set.tip_length, s.set.middle_length
    s.γ_l = pi/2 - width/2/radius
    s.γ_D = s.γ_l + width*(-2*tip_length + sqrt(2*middle_length^2 + 2*tip_length^2)) /
        (4*(middle_length - tip_length)) / radius
    s.kite_length_D = tip_length + (middle_length-tip_length) * (s.γ_D - s.γ_l) / (π/2 - s.γ_l)

    calc_inertia!(s)
    calc_pos_principal!(s)
    nothing
end

# include(joinpath(@__DIR__, "CreatePolars.jl"))
function KPSQ(kcu::KCU)
    set = kcu.set
    @assert set.foil_file != "" "No foil file specified in settings."
    open(joinpath(dirname(get_data_path()), set.foil_file), "r") do f
        lines = readlines(f)
        if !endswith(chomp(lines[1]), "polars created")
            error("No polars created for $(s.set.foil_file). Run scripts/create_polars.jl to create a polars file.")
        end
    end

    alphas, d_trailing_edge_angles, cl_matrix, cd_matrix, c_te_matrix = deserialize(joinpath(dirname(get_data_path()), set.polar_file))
    replace_nan!(cl_matrix)
    replace_nan!(cd_matrix)
    replace_nan!(c_te_matrix)
    cl_struct = extrapolate(scale(interpolate(cl_matrix, BSpline(Quadratic())), alphas, d_trailing_edge_angles), NaN)
    cd_struct = extrapolate(scale(interpolate(cd_matrix, BSpline(Quadratic())), alphas, d_trailing_edge_angles), NaN)
    c_te_struct = extrapolate(scale(interpolate(c_te_matrix, BSpline(Quadratic())), alphas, d_trailing_edge_angles), NaN)
    cl_interp(a, d) = cl_struct(a, d)
    cd_interp(a, d) = cd_struct(a, d)
    c_te_interp(a, d) = c_te_struct(a, d)
    
    if set.winch_model == "TorqueControlledMachine"
        s = KPSQ{SimFloat, Vector{SimFloat}, 3*(set.segments + 1)}(
            set=kcu.set, 
            motors=[TorqueControlledMachine(set) for _ in 1:3],
            cl_interp = cl_interp,
            cd_interp = cd_interp,
            c_te_interp = c_te_interp,)
        s.torque_control = true
    else
        s = KPSQ{SimFloat, Vector{SimFloat}, 3*(set.segments + 1)}(
            set=kcu.set, 
            motors=[AsyncMachine(set) for _ in 1:3],
            cl_interp = cl_interp,
            cd_interp = cd_interp,
            c_te_interp = c_te_interp,)
        s.torque_control = false
    end
    clear!(s)    
    return s
end

function calc_kite_ref_frame!(s::KPSQ, E, C, D)
    P_c = 0.5 .* (C+D)
    s.e_y .= normalize(C - D)
    s.e_z .= normalize(E - P_c)
    s.e_x .= cross(s.e_y, s.e_z)
    return nothing
end

function calc_tether_elevation(s::KPSQ)
    KiteUtils.calc_elevation(s.pos[6])
end

function calc_tether_azimuth(s::KPSQ)
    KiteUtils.azimuth_east(s.pos[6])
end

function update_sys_state!(ss::SysState, s::KPSQ, zoom=1.0)
    ss.time = s.t_0
    pos = cat(s.get_pos(), s.get_kite_pos(); dims=2)
    P = s.i_C + 1
    for i in 1:P
        ss.X[i] = pos[1, i] * zoom
        ss.Y[i] = pos[2, i] * zoom
        ss.Z[i] = pos[3, i] * zoom
    end
    # TODO
    # ss.kite_acc      .= kite_acc(s)
    # ss.left_tether_vel = tether_vel[1]
    # ss.right_tether_vel = tether_vel[2]
    ss.acc = norm(s.get_kite_acc())
    ss.orient .= quaternion_multiply(s.Q_p_b, s.get_Q_p_w())
    ss.elevation = calc_elevation(s)
    ss.azimuth = calc_azimuth(s)
    ss.force = tether_force(s)[3]
    ss.heading = calc_heading_y(s.get_e_x())
    ss.course = calc_course(s)
    ss.v_app = norm(s.v_apparent)
    ss.l_tether = s.tether_lengths[3]
    ss.v_reelout = s.get_tether_vel()[3]
    ss.depower = rad2deg(sum(s.get_trailing_edge_angle()))
    ss.steering = rad2deg(diff(s.get_trailing_edge_angle())[1])
    ss.vel_kite .= s.vel_kite
    nothing
end

function SysState(s::KPSQ, zoom=1.0) # TODO: add left and right lines, stop using getters and setters
    isnothing(s.integrator) && @warn "Initialize kite first!"
    pos = cat(s.get_pos(), s.get_kite_pos(); dims=2)
    P = s.i_C + 1
    X = zeros(MVector{P, MyFloat})
    Y = zeros(MVector{P, MyFloat})
    Z = zeros(MVector{P, MyFloat})
    for i in 1:P
        X[i] = pos[1, i] * zoom
        Y[i] = pos[2, i] * zoom
        Z[i] = pos[3, i] * zoom
    end
    
    orient = MVector{4, Float32}(calc_orient_quat(s))
    elevation = calc_elevation(s)
    azimuth = calc_azimuth(s)
    forces = s.get_tether_force()
    heading = calc_heading_y(s.get_e_x())
    course = calc_course(s)
    v_app_norm = norm(s.v_apparent)
    t_sim = 0
    depower = rad2deg(sum(s.get_trailing_edge_angle()))
    steering = rad2deg(diff(s.get_trailing_edge_angle())[1])
    ss = SysState{P}()
    ss.time = s.t_0
    ss.t_sim = t_sim
    ss.orient .= orient
    ss.elevation = elevation
    ss.azimuth = azimuth
    ss.l_tether = s.tether_lengths[3]
    ss.v_reelout = s.get_tether_vel()[3]
    ss.force = forces[3]
    ss.depower = depower
    ss.steering = steering
    ss.heading = heading
    ss.course = course
    ss.v_app = v_app_norm
    ss.vel_kite .= s.vel_kite
    ss.X = X
    ss.Y = Y
    ss.Z = Z
    ss
end


# function calc_heading(e_x, pos_kite)
#     # turn s.e_x by -azimuth around global z-axis and then by elevation around global y-axis
#     vec = rotate_around_y(rotate_around_z(e_x, -KiteUtils.azimuth_east(pos_kite)), -KiteUtils.calc_elevation(pos_kite))
#     heading = atan(-vec[2], vec[3])
#     return heading
# end
# @register_symbolic calc_heading(e_x, pos_kite)

function calc_heading_y(e_x)
    return atan(-e_x[2]/-e_x[1])
end
@register_symbolic calc_heading_y(e_x)

"""
Initialises the integrator of the model.

Parameters:
- s:     an instance of a 3 line kite model
- prn: if set to true, print the detailed solver results
- torque_control: wether or not to use torque control

Returns:
Nothing.
"""
function init_sim!(s::KPSQ; prn=false, torque_control=s.torque_control, 
        init_set_values=s.init_set_values, ϵ=s.ϵ, flap_damping=s.flap_damping, force_new_sys=false)
    dt = 1/s.set.sample_freq
    tspan   = (0.0, dt) 
    solver = TRBDF2( # https://docs.sciml.ai/SciMLBenchmarksOutput/stable/#Results
        autodiff=AutoFiniteDiff()
    ) # the best
    # solver = QBDF( # https://docs.sciml.ai/SciMLBenchmarksOutput/stable/#Results
    #     autodiff=AutoFiniteDiff()
    # )
    set_hash = struct_hash(s.set)
    measure_hash = struct_hash(s.measure)
    new_pos = s.last_measure_hash != measure_hash
    new_sys = force_new_sys ||
                    isnothing(s.prob) || 
                    s.torque_control != torque_control || 
                    s.last_set_hash != set_hash ||
                    s.ϵ != ϵ ||
                    s.flap_damping != flap_damping
    s.torque_control = torque_control
    s.ϵ = ϵ
    s.flap_damping = flap_damping

    if new_sys || true
        if prn; println("initializing with new model and new pos"); end
        clear!(s)
        init = true
        sys, u0map, p0map = model!(s; init)
        s.prob = ODEProblem(sys, u0map, tspan, p0map)
        s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
        generate_getters!(s; init)
        init_distance!(s)
    elseif new_pos
        if prn; println("initializing with last model and new pos"); end
        init_distance!(s)
    else
        if prn; println("initializing with last model and last pos"); end
        OrdinaryDiffEqCore.reinit!(s.integrator, s.prob.u0)
    end

    s.last_init_elevation = s.set.elevation
    s.last_init_tether_length = s.set.l_tether
    s.last_set_hash = set_hash
    s.last_measure_hash = measure_hash
    s.init_set_values .= init_set_values
    return nothing
end

function generate_getters!(s; init=false)
    sys = s.simple_sys

    if !init
        set_Q_p_w = setu(s.prob, sys.Q_p_w)
        set_ω_p = setu(s.prob, sys.ω_p)  
        set_kite_pos = setu(s.prob, sys.kite_pos)
        set_kite_vel = setu(s.prob, sys.kite_vel)
        set_pos = setu(s.prob, sys.pos[:, 4:s.i_A-1])
        set_vel = setu(s.prob, sys.vel[:, 4:s.i_A-1])
        set_trailing_edge_angle = setu(s.prob, sys.trailing_edge_angle)
        set_trailing_edge_ω = setu(s.prob, sys.trailing_edge_ω)
        set_gust_factor = setu(s.prob, sys.gust_factor)
        set_tether_length = setu(s.prob, sys.tether_length)
        set_tether_vel = setu(s.prob, sys.tether_vel)

        s.set_Q_p_w               = (prob, val) -> set_Q_p_w(prob, val)
        s.set_ω_p                 = (prob, val) -> set_ω_p(prob, val)  
        s.set_kite_pos            = (prob, val) -> set_kite_pos(prob, val)
        s.set_kite_vel            = (prob, val) -> set_kite_vel(prob, val)
        s.set_pos                 = (prob, val) -> set_pos(prob, val)
        s.set_vel                 = (prob, val) -> set_vel(prob, val)
        s.set_trailing_edge_angle = (prob, val) -> set_trailing_edge_angle(prob, val)
        s.set_trailing_edge_ω     = (prob, val) -> set_trailing_edge_ω(prob, val)
        s.set_gust_factor         = (prob, val) -> set_gust_factor(prob, val)
        s.set_tether_length       = (prob, val) -> set_tether_length(prob, val)
        s.set_tether_vel          = (prob, val) -> set_tether_vel(prob, val)
    end

    set_v_wind_gnd = setp(s.integrator, sys.v_wind_gnd)
    set_v_wind = setp(s.integrator, sys.v_wind)
    set_set_values = setu(s.integrator, sys.set_values)

    s.set_v_wind_gnd          = val -> set_v_wind_gnd(s.integrator, val)
    s.set_v_wind              = val -> set_v_wind(s.integrator, val)
    s.set_set_values          = val -> set_set_values(s.integrator, val)

    get_pos = getu(s.integrator, sys.pos)
    get_vel = getu(s.integrator, sys.vel)
    get_acc = getu(s.integrator, sys.acc)
    get_trailing_edge_angle = getu(s.integrator, sys.trailing_edge_angle)
    get_trailing_edge_α = getu(s.integrator, sys.trailing_edge_α)
    get_kite_pos = getu(s.integrator, sys.kite_pos)
    get_kite_vel = getu(s.integrator, sys.kite_vel)
    get_kite_acc = getu(s.integrator, sys.kite_acc)
    get_R_b_w = getu(s.integrator, sys.R_b_w)
    get_Q_p_w = getu(s.integrator, sys.Q_p_w)
    get_e_x = getu(s.integrator, sys.e_x)
    get_e_y = getu(s.integrator, sys.e_y)
    get_e_z = getu(s.integrator, sys.e_z)
    get_tether_force = getu(s.integrator, sys.tether_force)
    get_tether_length = getu(s.integrator, sys.tether_length)
    get_tether_vel = getu(s.integrator, sys.tether_vel)
    get_tether_acc = getu(s.integrator, sys.tether_acc)
    get_kite_force = getu(s.integrator, sys.total_kite_force)
    get_kite_torque_p = getu(s.integrator, sys.torque_p)
    get_heading = getu(s.integrator, sys.heading_y)
    get_ω_b = getu(s.integrator, sys.ω_b)
    get_α_b = getu(s.integrator, sys.α_b)
    get_force = getu(s.integrator, sys.force)
    get_e_te_A = getu(s.integrator, sys.e_te_A)
    get_e_te_B = getu(s.integrator, sys.e_te_B)
    get_elevation_vel = getu(s.integrator, sys.elevation_vel)
    get_elevation_acc = getu(s.integrator, sys.elevation_acc)
    get_azimuth_vel = getu(s.integrator, sys.azimuth_vel)
    get_azimuth_acc = getu(s.integrator, sys.azimuth_acc)
    get_distance_acc = getu(s.integrator, sys.distance_acc)

    s.get_pos                 = () -> get_pos(s.integrator)
    s.get_vel                 = () -> get_vel(s.integrator)
    s.get_acc                 = () -> get_acc(s.integrator)
    s.get_trailing_edge_angle = () -> get_trailing_edge_angle(s.integrator)
    s.get_trailing_edge_α     = () -> get_trailing_edge_α(s.integrator)
    s.get_kite_pos            = () -> get_kite_pos(s.integrator)
    s.get_kite_vel            = () -> get_kite_vel(s.integrator)
    s.get_kite_acc            = () -> get_kite_acc(s.integrator)
    s.get_R_b_w               = () -> get_R_b_w(s.integrator)
    s.get_Q_p_w               = () -> get_Q_p_w(s.integrator)
    s.get_e_x                 = () -> get_e_x(s.integrator)
    s.get_e_y                 = () -> get_e_y(s.integrator)
    s.get_e_z                 = () -> get_e_z(s.integrator)
    s.get_tether_force        = () -> get_tether_force(s.integrator)
    s.get_tether_length       = () -> get_tether_length(s.integrator)
    s.get_tether_vel          = () -> get_tether_vel(s.integrator)
    s.get_tether_acc          = () -> get_tether_acc(s.integrator)
    s.get_kite_force          = () -> get_kite_force(s.integrator)
    s.get_kite_torque_p       = () -> get_kite_torque_p(s.integrator)
    s.get_heading             = () -> get_heading(s.integrator)
    s.get_ω_b                 = () -> get_ω_b(s.integrator)
    s.get_α_b                 = () -> get_α_b(s.integrator)
    s.get_force               = () -> get_force(s.integrator)
    s.get_e_te_A              = () -> get_e_te_A(s.integrator)
    s.get_e_te_B              = () -> get_e_te_B(s.integrator)
    s.get_elevation_vel       = () -> get_elevation_vel(s.integrator)
    s.get_elevation_acc       = () -> get_elevation_acc(s.integrator)
    s.get_azimuth_vel         = () -> get_azimuth_vel(s.integrator)
    s.get_azimuth_acc         = () -> get_azimuth_acc(s.integrator)
    s.get_distance_acc        = () -> get_distance_acc(s.integrator)
    nothing
end

function next_step!(s::KPSQ; set_values=nothing, v_wind_gnd=s.set.v_wind, upwind_dir=-pi/2, dt=1/s.set.sample_freq)
    set_v_wind_ground!(s, calc_height(s), v_wind_gnd; upwind_dir)
    if !isnothing(set_values)
        s.set_set_values(set_values)
    end
    s.set_v_wind_gnd(s.v_wind_gnd)
    s.set_v_wind(s.v_wind)
    s.t_0 = s.integrator.t
    OrdinaryDiffEqCore.step!(s.integrator, dt, true)
    if !successful_retcode(s.integrator.sol)
        println("Return code for solution: ", s.integrator.sol.retcode)
    end
    @assert successful_retcode(s.integrator.sol)
    s.integrator.t
end

# function calc_pre_tension(s::KPSQ)
#     forces = spring_forces(s)
#     avg_force = 0.0
#     for i in 1:s.i_A
#         avg_force += forces[i]
#     end
#     avg_force /= s.i_A
#     res = avg_force/s.c_spring[3]
#     if res < 0.0 res = 0.0 end
#     if isnan(res) res = 0.0 end
#     return res + 1.0
# end

"""
    unstretched_length(s::KPSQ)

Getter for the unstretched tether reel-out length (at zero force).
"""
function unstretched_length(s::KPSQ) s.tether_lengths[3] end

"""
Only on x and z axis, tether start and end point laying on x axis
"""
function generate_tether!(pos, d, segments, tether_length, total_angle)
    segment_angle = total_angle / segments
    pos[:, 1] .= 0
    initial_angle = -total_angle / 2 + segment_angle/2
    for i in 2:segments+1
        angle = initial_angle + (i - 2) * segment_angle
        pos[1, i] = pos[1, i-1] + tether_length/segments * cos(angle)
        pos[2, i] = 0.0
        pos[3, i] = pos[3, i-1] + tether_length/segments * sin(angle)
    end
    dx = pos[1, end] - pos[1, 1]
    dz = pos[3, end] - pos[3, 1]
    d[] = sqrt(dx^2 + dz^2)
    nothing
end

"""
Rotate a 3d matrix by a quaternion rotation
"""
function rotate!(pos, quat)
    for i in eachindex(pos[1, :])
        pos[:, i] .= quat * pos[:, i]
    end
    nothing
end

"""
Generate a 2d tether given a certain distance and length.
"""
function tether_from_distance_length!(pos, distance::SimFloat, tether_length::SimFloat, segments, quat)
    d = Ref(0.0)
    cost = Ref(0.0)
    d[] = 0.0
    cost[] = 0.0
    function f_zero!(total_angle)
        generate_tether!(pos, d, segments, tether_length, total_angle)
        cost[] = d[] - distance
        return cost[]
    end
    Roots.find_zero(f_zero!, (0, 2π); atol=1e-6)
    rotate!(pos, quat)
    nothing
end

"""
Rotation from vector u to vector v
"""
function quaternion_rotation(u, v)
    d = dot(u, v)
    w = cross(u, v)
    return QuatRotation(d + sqrt(d * d + dot(w, w)), w...)
end

"""
Kite pos: C flap pos - D flap pos - middle tether pos \
Kite vel: C flap vel - D flap vel - middle tether vel \
Tether vel: left - middle - right tether vel

Return: an expected vel for all kite pos
"""
function calc_expected_pos_vel(s::KPSQ, kite_pos1, kite_pos2, kite_pos3, kite_vel, tether_vel, tether_length, tether_force, c_spring) # TODO: remove the 123 caused by this issue: https://github.com/SciML/ModelingToolkit.jl/issues/3003
    kite_pos = [kite_pos1, kite_pos2, kite_pos3]
    s.expected_tether_pos_vel_buffer .= 0.0
    expected_pos = @views s.expected_tether_pos_vel_buffer[1, :, :]
    expected_vel = @views s.expected_tether_pos_vel_buffer[2, :, :]
    J, y, x, segments = s.J_buffer, s.y_buffer, s.x_buffer, s.set.segments
    distance = norm(kite_pos)
    stretched_tether_length = tether_length + tether_force / (c_spring/tether_length)

    if any(isa.((kite_pos1, kite_pos2, kite_pos3, kite_vel, tether_vel, tether_length, tether_force, c_spring), ForwardDiff.Dual)) ||
            !all(0.0 .< (distance, stretched_tether_length) .< Inf) || 
            distance >= stretched_tether_length ||
            distance <= 0.1stretched_tether_length
        expected_pos .= NaN
        expected_vel .= NaN
        return s.expected_tether_pos_vel_buffer
    end

    quat = quaternion_rotation([1, 0, 0], kite_pos)
    function f_jac!(dx, x)
        tether_from_distance_length!(expected_pos, x[1], x[2], segments, quat)
        dx .= vec(expected_pos)
        nothing
    end

    x[1] = distance
    x[2] = stretched_tether_length

    backend = ADTypes.AutoFiniteDiff()
    if isnothing(s.prep)
        s.prep = prepare_jacobian(f_jac!, y, backend, x)
    end
    DifferentiationInterface.jacobian!(f_jac!, y, J, s.prep, backend, x)
    expected_vel .= reshape(J * [kite_vel, tether_vel], size(expected_vel))
    return s.expected_tether_pos_vel_buffer
end
const FD = ForwardDiff.Dual
function calc_expected_pos_vel(s::KPSQ, _::FD, _::FD, _::FD, _::FD, _::FD, _::FD, _::FD, _::FD) # dummy function for forwarddiff compatibility
    s.expected_tether_pos_vel_buffer .= NaN
    return s.expected_tether_pos_vel_buffer
end
@register_array_symbolic calc_expected_pos_vel(s::KPSQ, kite_pos1, kite_pos2, kite_pos3, kite_vel, tether_vel, tether_length, tether_force, c_spring) begin
    size = size(s.expected_tether_pos_vel_buffer)
    eltype = SimFloat
end


# =================== getter functions ====================================================

"""
    calc_height(s::KPSQ)

Determine the height of the topmost kite particle above ground.
"""
function calc_height(s::KPSQ)
    pos_kite(s)[3]
end

"""
    pos_kite(s::KPSQ)

Return the position of the kite (top particle).
"""
function pos_kite(s::KPSQ)
    s.kite_pos
end

"""
    kite_ref_frame(s::KPSQ; one_point=false)

Returns a tuple of the x, y, and z vectors of the kite reference frame.
The parameter one_point is not used in this model.
"""
function kite_ref_frame(s::KPSQ; one_point=false)
    s.get_e_x(), s.get_e_y(), s.get_e_z()
end

"""
    tether_force(s::KPSQ)

Return the absolute value of the force at the winch as calculated during the last timestep. 
"""
function tether_force(s::KPSQ) s.get_tether_force() end

# ====================== helper functions ====================================

function struct_hash(st)
    fields = fieldnames(typeof(st))
    h = zero(UInt)
    for field in fields
        field_value = getfield(st, field)
        h = hash(field_value, h)
    end
    return h
end

function replace_nan!(matrix)
    rows, cols = size(matrix)
    distance = max(rows, cols)
    for i in 1:rows
        for j in 1:cols
            if isnan(matrix[i, j])
                neighbors = []
                for d in 1:distance
                    found = false
                    if i-d >= 1 && !isnan(matrix[i-d, j]);
                        push!(neighbors, matrix[i-d, j])
                        found = true
                    end
                    if i+d <= rows && !isnan(matrix[i+d, j])
                        push!(neighbors, matrix[i+d, j])
                        found = true
                    end
                    if j-d >= 1 && !isnan(matrix[i, j-d])
                        push!(neighbors, matrix[i, j-d])
                        found = true
                    end
                    if j+d <= cols && !isnan(matrix[i, j+d])
                        push!(neighbors, matrix[i, j+d])
                        found = true
                    end
                    if found; break; end
                end
                if !isempty(neighbors)
                    matrix[i, j] = sum(neighbors) / length(neighbors)
                end
            end
        end
    end
    return nothing
end
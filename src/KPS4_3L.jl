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
    distance::MeasureFloat      = 0.0
    distance_acc::MeasureFloat  = 0.0
    winch_torque::MVector{3, MeasureFloat}  = zeros(3)
    tether_length::MVector{3, MeasureFloat} = zeros(3)
    tether_vel::MVector{3, MeasureFloat}    = zeros(3)
    tether_acc::MVector{3, MeasureFloat}    = zeros(3)
    azimuth_left::MeasureFloat       = 0.0
    d_azimuth_left::MeasureFloat     = 0.0
    dd_azimuth_left::MeasureFloat    = 0.0
    azimuth_right::MeasureFloat       = 0.0
    d_azimuth_right::MeasureFloat     = 0.0
    dd_azimuth_right::MeasureFloat    = 0.0
    elevation_left::MeasureFloat     = 0.0
    d_elevation_left::MeasureFloat   = 0.0
    dd_elevation_left::MeasureFloat  = 0.0
    elevation_right::MeasureFloat     = 0.0
    d_elevation_right::MeasureFloat   = 0.0
    dd_elevation_right::MeasureFloat  = 0.0
end

const Getter = Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing}

"""
    mutable struct KPS4_3L{S, T, P, Q, SP} <: AbstractKiteModel

State of the kite power system, using a 3 point kite model and three steering lines to the ground. Parameters:
- S: Scalar type, e.g. SimFloat
  In the documentation mentioned as Any, but when used in this module it is always SimFloat and not Any.
- V: Vector type, e.g. KVec3
- P: number of tether points of the system, 3segments+3
Normally a user of this package will not have to access any of the members of this type directly,
use the input and output functions instead.

$(TYPEDFIELDS)
"""
@with_kw mutable struct KPS4_3L{S, V, P} <: AbstractKiteModel # TODO: subdivide in kite-changing fields and non-kite-changing fields. combine this with fast caching
    "Reference to the settings struct"
    set::Settings
    "Reference to the settings hash"
    set_hash::UInt64        = 0
    "The last initial elevation"
    last_init_elevation::S     = 0.0
    "The last initial tether length"
    last_init_tether_length::S = 0.0
    "Reference to the last settings hash"
    last_set_hash::UInt64   = 0
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
    pos::Matrix{S} = zeros(3, P)
    "unstressed segment lengths of the three tethers [m]"
    segment_lengths::V =           zeros(S, 3)
    "relative start time of the current time interval"
    t_0::S =               0.0
    "unstretched tether length"
    tether_lengths::V =          zeros(S, 3)
    "air density at the height of the kite"
    rho::S =               0.0
    "multiplier for the damping of all movement"
    damping_coeff::S =  50.0
    "tether masses"
    masses::V         = zeros(P)
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
    α_l::S =     0.0
    "Angle of point C"
    α_D::S =     0.0
    "Kite length at point C"
    kite_length_D::S =     0.0
    "Simplified system of the mtk model"
    simple_sys::Union{ModelingToolkit.ODESystem, Nothing} = nothing
    "Velocity of the kite"
    vel_kite::V =           zeros(S, 3)
    "Initial torque or speed set values"
    init_set_values::V =    zeros(S, 3)
    "Smooth sign constant"
    ϵ::S =      1e-3
    "Relative damping for flaps"
    flap_damping::S     = 0.75
    "Measured data points used to create an initial state"
    measure::Measurements = Measurements()
    "Buffer for expected kite pos and vel, defined as (pos, vel), (x, y, z), (value)"
    expected_tether_pos_vel_buffer::Array{S, 3} = zeros(2, 3, (P ÷ 3))
    "Buffers for Jacobian for expected tether point velocities relative to winch velocities and kite velocities"
    J_buffer::Matrix{S} = zeros(P, 2)
    "Makes autodiff faster"
    prep::Union{Nothing, Any} = nothing
    "Buffer for jacobian y values (vectorized velocities)"
    y_buffer::V = zeros(P)
    "Buffer for jacobian x values (angle, distance)"
    x_buffer::V = zeros(2)
    "Inertia around kite x y and z axis"
    I_kite::V = zeros(3)
    "Damping of the kite rotation"
    orient_damping::S = zero(S)
    "rotation from kite body frame to kite principal frame"
    R_b_p::Matrix{S} = zeros(3, 3)
    "translation from kite point to circle center in body frame along positive z"
    circle_center_t::V = zeros(3)
    "center of mass of the kite"
    kite_pos::V = zeros(3)
    "quaternion orientation of the kite"
    q::V = zeros(4)
    "aero points center of mass positions in principal frame"
    seg_com_pos_p::Matrix{S} = zeros(3, 2set.aero_surfaces)
    "aero points center of pressure positions in principal frame"
    seg_cop_pos_p::Matrix{S} = zeros(3, 2set.aero_surfaces)
    "cop pos in body frame"
    seg_cop_pos_b::Matrix{S} = zeros(3, 2set.aero_surfaces)
    "last left tether point position in principal frame"
    pos_A_p::V = zeros(3)
    "last right tether point position in principal frame"
    pos_B_p::V = zeros(3)
    "last middle tether point position in principal frame"
    pos_C_p::V = zeros(3)
    "rotation from kite body frame to world frame"
    R_b_w::Matrix{S} = zeros(3, 3)
    "translation from kite_point to point C along positive body z frame"
    C_t::V = zeros(3)
    "mass of each kite segment"
    seg_mass::V = zeros(2set.aero_surfaces)
    "A point projected onto kite in z-axis in body frame"
    pos_D_b::V = zeros(3)
    "B point projected onto kite in z-axis in body frame"
    pos_E_b::V = zeros(3)

    # set_values_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    v_wind_gnd_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    v_wind_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    prob::Union{OrdinaryDiffEqCore.ODEProblem, Nothing} = nothing
    set_set_values::Union{SymbolicIndexingInterface.MultipleSetters, Nothing} = nothing
    get_pos::Getter = nothing
    get_kite_pos::Getter = nothing
    get_trailing_edge_angle::Getter = nothing
    get_trailing_edge_α::Getter = nothing
    get_kite_vel::Getter = nothing
    get_kite_acc::Getter = nothing
    get_q::Getter = nothing
    get_tether_forces::Getter = nothing
    get_tether_lengths::Getter = nothing
    get_tether_vels::Getter = nothing
    get_kite_force::Getter = nothing
    get_kite_torque_p::Getter = nothing
    get_heading::Getter = nothing
    integrator::Union{OrdinaryDiffEqCore.ODEIntegrator, Nothing} = nothing
    u0::V = [0.0]
end

function get_kite_torque_b(s)
    return s.R_b_p' * s.get_kite_torque_p(s.integrator)
end

"""
    clear!(s::KPS4_3L)

Initialize the kite power model.
"""
function clear!(s::KPS4_3L)
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
    s.α_l = π/2 - s.set.min_steering_line_distance/(2*s.set.radius)
    s.segment_lengths .= s.tether_lengths ./ s.set.segments
    s.i_A = s.set.segments*3+1
    s.i_B = s.set.segments*3+2
    s.i_C = s.set.segments*3+3
    s.rho = s.set.rho_0
    c_spring = s.set.e_tether * (s.set.d_tether/2000.0)^2 * pi
    s.c_spring .= [c_spring, c_spring, 2c_spring]
    s.damping .= (s.set.damping / s.set.c_spring) * s.c_spring
    init_masses!(s)
end

# include(joinpath(@__DIR__, "CreatePolars.jl"))
function KPS4_3L(kcu::KCU)
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
        s = KPS4_3L{SimFloat, Vector{SimFloat}, 3*(set.segments + 1)}(
            set=kcu.set, 
            motors=[TorqueControlledMachine(set) for _ in 1:3],
            cl_interp = cl_interp,
            cd_interp = cd_interp,
            c_te_interp = c_te_interp,)
        s.torque_control = true
    else
        s = KPS4_3L{SimFloat, Vector{SimFloat}, 3*(set.segments + 1)}(
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

function calc_kite_ref_frame!(s::KPS4_3L, E, C, D)
    P_c = 0.5 .* (C+D)
    s.e_y .= normalize(C - D)
    s.e_z .= normalize(E - P_c)
    s.e_x .= cross(s.e_y, s.e_z)
    return nothing
end

function calc_tether_elevation(s::KPS4_3L)
    KiteUtils.calc_elevation(s.pos[6])
end

function calc_tether_azimuth(s::KPS4_3L)
    KiteUtils.azimuth_east(s.pos[6])
end

function update_sys_state!(ss::SysState, s::KPS4_3L, zoom=1.0)
    ss.time = s.t_0
    pos = cat(s.get_pos(s.integrator), s.get_kite_pos(s.integrator); dims=2)
    P = s.i_C + 1
    for i in 1:P
        ss.X[i] = pos[1, i] * zoom
        ss.Y[i] = pos[2, i] * zoom
        ss.Z[i] = pos[3, i] * zoom
    end
    tether_vels = s.get_tether_vels(s.integrator)
    # TODO
    # ss.kite_acc      .= s.get_kite_acc(s.integrator)
    # ss.left_tether_vel = tether_vels[1]
    # ss.right_tether_vel = tether_vels[2]
    ss.acc = norm(s.get_kite_acc(s.integrator))
    ss.orient .= s.get_q(s.integrator)
    ss.elevation = calc_elevation(s)
    ss.azimuth = calc_azimuth(s)
    ss.force = tether_force(s)[3]
    ss.heading = calc_heading_y(s.e_x)
    ss.course = calc_course(s)
    ss.v_app = norm(s.v_apparent)
    ss.l_tether = s.tether_lengths[3]
    ss.v_reelout = s.get_tether_vels(s.integrator)[3]
    ss.depower = rad2deg(s.get_trailing_edge_angle(s.integrator)[1] + s.get_trailing_edge_angle(s.integrator)[2])
    ss.steering = rad2deg(s.get_trailing_edge_angle(s.integrator)[2] - s.get_trailing_edge_angle(s.integrator)[1])
    ss.vel_kite .= s.vel_kite
    nothing
end

function SysState(s::KPS4_3L, zoom=1.0) # TODO: add left and right lines, stop using getters and setters
    isnothing(s.integrator) && @warn "Initialize kite first!"
    generate_getters!(s)
    pos = cat(s.get_pos(s.integrator), s.get_kite_pos(s.integrator); dims=2)
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
    forces = tether_force(s)
    heading = calc_heading_y(s.e_x)
    course = calc_course(s)
    v_app_norm = norm(s.v_apparent)
    t_sim = 0
    depower = rad2deg(s.get_trailing_edge_angle(s.integrator)[1] + s.get_trailing_edge_angle(s.integrator)[2])
    steering = rad2deg(s.get_trailing_edge_angle(s.integrator)[2] - s.get_trailing_edge_angle(s.integrator)[1])
    ss = SysState{P}()
    ss.time = s.t_0
    ss.t_sim = t_sim
    ss.orient .= orient
    ss.elevation = elevation
    ss.azimuth = azimuth
    ss.l_tether = s.tether_lengths[3]
    ss.v_reelout = s.get_tether_vels(s.integrator)[3]
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
    init_sim!(s::KPS4_3L; damping_coeff=50.0, prn=false, torque_control=true)

Initialises the integrator of the model.

Parameters:
- s:     an instance of a 3 line kite model
- damping_coeff: amount of damping in the first steps
- prn: if set to true, print the detailed solver results
- torque_control: wether or not to use torque control

Returns:
Nothing.
"""
function init_sim!(s::KPS4_3L; damping_coeff=s.damping_coeff, prn=false, 
                   torque_control=s.torque_control, init_set_values=s.init_set_values, ϵ=s.ϵ, flap_damping=s.flap_damping)
    clear!(s)
    
    dt = 1/s.set.sample_freq
    tspan   = (0.0, dt) 
    solver = QNDF(autodiff=false) # https://docs.sciml.ai/SciMLBenchmarksOutput/stable/#Results
    new_inital_conditions =     (s.last_init_elevation != s.set.elevation || 
                                s.last_init_tether_length != s.set.l_tether) || 
                                !all(s.init_set_values .== init_set_values)
    s.set_hash = settings_hash(s.set)
    init_new_model = isnothing(s.prob) || 
                    s.torque_control != torque_control || 
                    s.last_set_hash != s.set_hash ||
                    s.damping_coeff != damping_coeff ||
                    s.ϵ != ϵ ||
                    s.flap_damping != flap_damping
    s.torque_control = torque_control
    s.damping_coeff = damping_coeff
    s.ϵ = ϵ
    s.flap_damping = flap_damping
    init_new_pos = new_inital_conditions && !isnothing(s.get_pos)

    dt0 = 1.0
    function stabilize_kite(s::KPS4_3L)
        set_values = copy(init_set_values)
        for _ in 0:dt:dt0
            next_step!(s; set_values, dt) # step to get stable state
            if (torque_control) set_values .= -tether_force(s) .* s.set.drum_radius end
        end
    end
    if init_new_model
        if prn; println("initializing with new model and new pos"); end
        pos, vel = init_pos_vel(s)
        sys, inputs = create_sys!(s, pos, vel)
        (s.simple_sys, _) = structural_simplify(sys, (inputs, []); simplify=true)
        s.prob = ODEProblem(s.simple_sys, nothing, tspan; fully_determined=true)
        s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
        stabilize_kite(s)
        s.u0 = deepcopy(s.integrator.u)
        OrdinaryDiffEqCore.reinit!(s.integrator, s.u0; t0=dt0, tf=dt0+dt)
    elseif init_new_pos
        if prn; println("initializing with last model and new pos"); end
        pos, vel = init_pos_vel(s)
        pos, vel = convert_pos_vel(s, pos, vel)
        defaults = vcat(
                    vcat([s.simple_sys.pos[j, i] => pos[j, i] for i in 1:s.i_A-1 for j in 1:3]), 
                    vcat([s.simple_sys.pos[j, i] => pos[j, i] for i in s.i_B+1:s.i_A for j in 1:3]),
                    vcat([s.simple_sys.tether_length[i] => s.tether_lengths[i] for i in 1:3]),
                        )
        s.prob = ODEProblem(s.simple_sys, defaults, tspan)
        OrdinaryDiffEqCore.reinit!(s.integrator, s.prob.u0)
        stabilize_kite(s)
        s.u0 = deepcopy(s.integrator.u)
        OrdinaryDiffEqCore.reinit!(s.integrator, s.u0; t0=dt0, tf=dt0+dt)
    else
        if prn; println("initializing with last model and last pos"); end
        OrdinaryDiffEqCore.reinit!(s.integrator, s.u0; t0=dt0, tf=dt0+dt)
        # s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
    end

    s.last_init_elevation = s.set.elevation
    s.last_init_tether_length = s.set.l_tether
    s.last_set_hash = s.set_hash
    s.init_set_values .= init_set_values
    return nothing
end

function generate_getters!(s)
    if isnothing(s.get_pos)
        s.v_wind_gnd_idx = parameter_index(s.integrator.f, :v_wind_gnd)
        s.v_wind_idx = parameter_index(s.integrator.f, :v_wind)
        s.set_set_values = setu(s.integrator.sol, s.simple_sys.set_values)
        s.get_pos = getu(s.integrator.sol, s.simple_sys.pos[:,:])
        s.get_trailing_edge_angle = getu(s.integrator.sol, s.simple_sys.trailing_edge_angle)
        s.get_trailing_edge_α = getu(s.integrator.sol, s.simple_sys.trailing_edge_α)
        s.get_kite_pos = getu(s.integrator.sol, s.simple_sys.kite_pos)
        s.get_kite_vel = getu(s.integrator.sol, s.simple_sys.kite_vel)
        s.get_kite_acc = getu(s.integrator.sol, s.simple_sys.kite_acc)
        s.get_q = getu(s.integrator.sol, s.simple_sys.q)
        s.get_tether_forces = getu(s.integrator.sol, s.simple_sys.tether_force)
        s.get_tether_lengths = getu(s.integrator.sol, s.simple_sys.tether_length)
        s.get_tether_vels = getu(s.integrator.sol, s.simple_sys.tether_vel)
        s.get_kite_force = getu(s.integrator.sol, s.simple_sys.total_kite_force)
        s.get_kite_torque_p = getu(s.integrator.sol, s.simple_sys.torque_p)
        s.get_heading = getu(s.integrator.sol, s.simple_sys.heading_y)
    end
    nothing
end

function next_step!(s::KPS4_3L; set_values=zeros(KVec3), v_wind_gnd=s.set.v_wind, upwind_dir=-pi/2, dt=1/s.set.sample_freq)
    set_v_wind_ground!(s, calc_height(s), v_wind_gnd; upwind_dir)
    generate_getters!(s)
    s.set_set_values(s.integrator, set_values)
    s.integrator.p[s.v_wind_gnd_idx] .= s.v_wind_gnd
    s.integrator.p[s.v_wind_idx] .= s.v_wind
    s.t_0 = s.integrator.t
    OrdinaryDiffEqCore.step!(s.integrator, dt, true)
    if !successful_retcode(s.integrator.sol)
        println("Return code for solution: ", s.integrator.sol.retcode)
    end
    @assert successful_retcode(s.integrator.sol)
    s.integrator.t
end

# function calc_pre_tension(s::KPS4_3L)
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
    unstretched_length(s::KPS4_3L)

Getter for the unstretched tether reel-out length (at zero force).
"""
function unstretched_length(s::KPS4_3L) s.tether_lengths[3] end

"""
    tether_length(s::KPS4_3L)

Calculate and return the real, stretched tether length.
"""
function tether_length(s::KPS4_3L)
    length = 0.0
    for i in 3:3:s.i_A-1
        length += norm(s.pos[i+3] - s.pos[i])
    end
    return length
end

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
    try
        Roots.find_zero(f_zero!, (0, 2π); atol=1e-6)
    catch e
        println("Distance: ", distance, " Tether length: ", tether_length)
    end
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
function calc_expected_pos_vel(s::KPS4_3L, kite_pos1, kite_pos2, kite_pos3, kite_vel, tether_vel, tether_length, tether_force, c_spring) # TODO: remove the 123 caused by this issue: https://github.com/SciML/ModelingToolkit.jl/issues/3003
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
function calc_expected_pos_vel(s::KPS4_3L, _::FD, _::FD, _::FD, _::FD, _::FD, _::FD, _::FD, _::FD) # dummy function for forwarddiff compatibility
    s.expected_tether_pos_vel_buffer .= NaN
    return s.expected_tether_pos_vel_buffer
end
@register_array_symbolic calc_expected_pos_vel(s::KPS4_3L, kite_pos1, kite_pos2, kite_pos3, kite_vel, tether_vel, tether_length, tether_force, c_spring) begin
    size = size(s.expected_tether_pos_vel_buffer)
    eltype = SimFloat
end


# =================== getter functions ====================================================

"""
    calc_height(s::KPS4_3L)

Determine the height of the topmost kite particle above ground.
"""
function calc_height(s::KPS4_3L)
    pos_kite(s)[3]
end

"""
    pos_kite(s::KPS4_3L)

Return the position of the kite (top particle).
"""
function pos_kite(s::KPS4_3L)
    s.kite_pos
end

"""
    kite_ref_frame(s::KPS4_3L; one_point=false)

Returns a tuple of the x, y, and z vectors of the kite reference frame.
The parameter one_point is not used in this model.
"""
function kite_ref_frame(s::KPS4_3L; one_point=false)
    s.e_x, s.e_y, s.e_z
end

"""
    tether_force(s::KPS4_3L)

Return the absolute value of the force at the winch as calculated during the last timestep. 
"""
function tether_force(s::KPS4_3L) s.get_tether_forces(s.integrator) end


# ====================== helper functions ====================================

function settings_hash(st)
    fields = fieldnames(typeof(st))
    fields = filter(x -> x != :l_tether && x != :elevation, fields)
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
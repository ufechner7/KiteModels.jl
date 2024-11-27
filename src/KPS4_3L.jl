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
configured in the file data/settings.yaml). The kite is modelled using 4 point masses and 3*n aerodynamic 
surfaces. The spring constant and the damping decrease with the segment length. The aerodynamic kite forces
are acting on three of the four kite point masses. 

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
    heading::MeasureFloat       = 0.0
    turn_rate::MeasureFloat     = 0.0
    turn_acc::MeasureFloat      = 0.0
    winch_torque::MVector{3, MeasureFloat}  = zeros(3)
    tether_length::MVector{3, MeasureFloat} = zeros(3)
    tether_vel::MVector{3, MeasureFloat}    = zeros(3)
    tether_acc::MVector{3, MeasureFloat}    = zeros(3)
    azimuth::MeasureFloat       = 0.0
    d_azimuth::MeasureFloat     = 0.0
    dd_azimuth::MeasureFloat    = 0.0
    elevation::MeasureFloat     = 0.0
    d_elevation::MeasureFloat   = 0.0
    dd_elevation::MeasureFloat  = 0.0
end

"""
    mutable struct KPS4_3L{S, T, P, Q, SP} <: AbstractKiteModel

State of the kite power system, using a 3 point kite model and three steering lines to the ground. Parameters:
- S: Scalar type, e.g. SimFloat
  In the documentation mentioned as Any, but when used in this module it is always SimFloat and not Any.
- T: Vector type, e.g. KVec3
- P: number of points of the system, segments+3
- Q: number of springs in the system, P-1
- SP: struct type, describing a spring
Normally a user of this package will not have to access any of the members of this type directly,
use the input and output functions instead.

$(TYPEDFIELDS)
"""
@with_kw mutable struct KPS4_3L{S, T, P, Q, SP} <: AbstractKiteModel
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
    "Iterations, number of calls to the function residual!"
    iter:: Int64 = 0
    "wind vector at the height of the kite" 
    v_wind::T =           zeros(S, 3)
    "wind vector at reference height" 
    v_wind_gnd::T =       zeros(S, 3)
    "wind vector used for the calculation of the tether drag"
    v_wind_tether::T =    zeros(S, 3)
    "apparent wind vector at the kite"
    v_apparent::T =       zeros(S, 3)
    "a copy of the residual one (pos,vel) for debugging and unit tests"    
    res1::SVector{P, T} = zeros(SVector{P, T})
    "a copy of the residual two (vel,acc) for debugging and unit tests"
    res2::SVector{P, T} = zeros(SVector{P, T})
    "a copy of the actual positions as output for the user"
    pos::SVector{P, T} = zeros(SVector{P, T})
    vel::SVector{P, T} = zeros(SVector{P, T})
    "unstressed segment lengths of the three tethers [m]"
    segment_lengths::T =           zeros(S, 3)
    "azimuth angle in radian; inital value is zero"
    psi::S =              zero(S)
    "relative start time of the current time interval"
    t_0::S =               0.0
    "unstretched tether length"
    tether_lengths::T =          zeros(S, 3)
    "air density at the height of the kite"
    rho::S =               0.0
    "multiplier for the damping of all movement"
    damping_coeff::S =  50.0
    "current masses, depending on the total tether length"
    masses::MVector{P, S}         = zeros(P)
    "vector of the springs, defined as struct"
    springs::MVector{Q, SP}       = zeros(SP, Q)
    "unit spring coefficient"
    c_spring::S = zero(S)
    "unit damping coefficient"
    damping::S = zero(S)
    "whether or not to use torque control instead of speed control"
    torque_control::Bool = false
    "x vector of kite reference frame"
    e_x::T =                 zeros(S, 3)
    "y vector of kite reference frame"
    e_y::T =                 zeros(S, 3)
    "z vector of kite reference frame"
    e_z::T =                 zeros(S, 3)
    "Point number of C flap connection point"
    num_flap_C::Int64 =           0
    "Point number of D flap connection point"
    num_flap_D::Int64 =           0
    "Point number of E"
    num_E::Int64 =           0
    "Point number of C"
    num_C::Int64 =           0
    "Point number of D"
    num_D::Int64 =           0
    "Point number of A"
    num_A::Int64 =           0
    "Angle of left tip"
    α_l::S =     0.0
    "Angle of right tip"
    α_r::S =     0.0
    "Angle of point C"
    α_C::S =     0.0
    "Kite length at point C"
    kite_length_C::S =     0.0
    "Solution of the steady state problem"
    steady_sol::Union{SciMLBase.NonlinearSolution, Nothing} = nothing
    "Simplified system of the mtk model"
    simple_sys::Union{ModelingToolkit.ODESystem, Nothing} = nothing
    "Velocity of the kite"
    vel_kite::T =           zeros(S, 3)
    "Initial torque or speed set values"
    init_set_values::T =    zeros(S, 3)
    "Smooth sign constant"
    ϵ::S =      1e-3
    "Relative damping for flaps"
    flap_damping::S     = 0.75
    "Measured data points used to create an initial state"
    measure::Measurements = Measurements()

    # set_values_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    v_wind_gnd_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    v_wind_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    prob::Union{OrdinaryDiffEqCore.ODEProblem, Nothing} = nothing
    set_set_values::Union{SymbolicIndexingInterface.MultipleSetters, Nothing} = nothing
    get_pos::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_flap_angle::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_flap_acc::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_vel_kite::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_winch_forces::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_tether_lengths::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_tether_vels::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_L_C::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_L_D::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_D_C::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_D_D::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_heading::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    integrator::Union{Sundials.CVODEIntegrator, OrdinaryDiffEqCore.ODEIntegrator, Nothing} = nothing
    u0:: Vector{SimFloat} = [0.0]
end

"""
    clear!(s::KPS4_3L)

Initialize the kite power model.
"""
function clear!(s::KPS4_3L)
    s.iter = 0
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
    s.α_r = π/2 + s.set.min_steering_line_distance/(2*s.set.radius)
    s.segment_lengths .= s.tether_lengths ./ s.set.segments
    s.num_flap_C = s.set.segments*3+3-2
    s.num_flap_D = s.set.segments*3+3-1
    s.num_E = s.set.segments*3+3
    s.num_C = s.set.segments*3+3+1
    s.num_D = s.set.segments*3+3+2
    s.num_A = s.set.segments*3+3+3
    s.rho = s.set.rho_0
    s.c_spring = s.set.e_tether * (s.set.d_tether/2000.0)^2 * pi
    s.damping = (s.set.damping / s.set.c_spring) * s.c_spring
    init_masses!(s)
    init_springs!(s)
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

    alphas, d_flap_angles, cl_matrix, cd_matrix, c_te_matrix = deserialize(joinpath(dirname(get_data_path()), set.polar_file))
    replace_nan!(cl_matrix)
    replace_nan!(cd_matrix)
    replace_nan!(c_te_matrix)
    cl_struct = extrapolate(scale(interpolate(cl_matrix, BSpline(Quadratic())), alphas, d_flap_angles), NaN)
    cd_struct = extrapolate(scale(interpolate(cd_matrix, BSpline(Quadratic())), alphas, d_flap_angles), NaN)
    c_te_struct = extrapolate(scale(interpolate(c_te_matrix, BSpline(Quadratic())), alphas, d_flap_angles), NaN)
    cl_interp(a, d) = cl_struct(a, d)
    cd_interp(a, d) = cd_struct(a, d)
    c_te_interp(a, d) = c_te_struct(a, d)
    
    if set.winch_model == "TorqueControlledMachine"
        s = KPS4_3L{SimFloat, KVec3, set.segments*3+2+KITE_PARTICLES_3L, set.segments*3+KITE_SPRINGS_3L, SP}(
            set=kcu.set, 
            motors=[TorqueControlledMachine(set) for _ in 1:3],
            cl_interp = cl_interp,
            cd_interp = cd_interp,
            c_te_interp = c_te_interp,)
        s.torque_control = true
    else
        s = KPS4_3L{SimFloat, KVec3, set.segments*3+2+KITE_PARTICLES_3L, set.segments*3+KITE_SPRINGS_3L, SP}(
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
    pos = s.pos
    P = length(pos)
    for i in 1:P
        ss.X[i] = pos[i][1] * zoom
        ss.Y[i] = pos[i][2] * zoom
        ss.Z[i] = pos[i][3] * zoom
    end
    ss.orient .= calc_orient_quat(s)
    ss.elevation = calc_elevation(s)
    ss.azimuth = calc_azimuth(s)
    ss.force = winch_force(s)[3]
    ss.heading = calc_heading(s.e_x, s.pos[s.num_E])
    ss.course = calc_course(s)
    ss.v_app = norm(s.v_apparent)
    ss.l_tether = s.tether_lengths[3]
    ss.v_reelout = s.get_tether_vels(s.integrator)[3]
    ss.depower = rad2deg(s.get_flap_angle(s.integrator)[1] + s.get_flap_angle(s.integrator)[2])
    ss.steering = rad2deg(s.get_flap_angle(s.integrator)[2] - s.get_flap_angle(s.integrator)[1])
    ss.vel_kite .= s.vel_kite
    nothing
end

function SysState(s::KPS4_3L, zoom=1.0) # TODO: add left and right lines, stop using getters and setters
    pos = s.pos
    P = length(pos)
    X = zeros(MVector{P, MyFloat})
    Y = zeros(MVector{P, MyFloat})
    Z = zeros(MVector{P, MyFloat})
    for i in 1:P
        X[i] = pos[i][1] * zoom
        Y[i] = pos[i][2] * zoom
        Z[i] = pos[i][3] * zoom
    end
    
    orient = MVector{4, Float32}(calc_orient_quat(s))
    elevation = calc_elevation(s)
    azimuth = calc_azimuth(s)
    forces = winch_force(s)
    heading = calc_heading(s.e_x, s.pos[s.num_E])
    course = calc_course(s)
    v_app_norm = norm(s.v_apparent)
    t_sim = 0
    depower = rad2deg(s.get_flap_angle(s.integrator)[1] + s.get_flap_angle(s.integrator)[2])
    steering = rad2deg(s.get_flap_angle(s.integrator)[2] - s.get_flap_angle(s.integrator)[1])
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


function calc_heading(e_x, pos_kite)
    # turn s.e_x by -azimuth around global z-axis and then by elevation around global y-axis
    vec = rotate_in_xz(rotate_in_yx(e_x, -KiteUtils.azimuth_east(pos_kite)), -KiteUtils.calc_elevation(pos_kite))
    heading = atan(-vec[2], vec[3])
    return heading
end
@register_symbolic calc_heading(e_x, pos_kite)

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
            if (torque_control) set_values .= -winch_force(s) .* s.set.drum_radius end
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
                    vcat([s.simple_sys.pos[j, i] => pos[j, i] for i in 1:s.num_flap_C-1 for j in 1:3]), 
                    vcat([s.simple_sys.pos[j, i] => pos[j, i] for i in s.num_flap_D+1:s.num_A for j in 1:3]),
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
    update_state!(s)
    return nothing
end


function next_step!(s::KPS4_3L; set_values=zeros(KVec3), v_wind_gnd=s.set.v_wind, upwind_dir=-pi/2, dt=1/s.set.sample_freq)
    s.iter = 0
    set_v_wind_ground!(s, calc_height(s), v_wind_gnd; upwind_dir)
    if isnothing(s.get_pos)
        s.v_wind_gnd_idx = parameter_index(s.integrator.f, :v_wind_gnd)
        s.v_wind_idx = parameter_index(s.integrator.f, :v_wind)
        s.set_set_values = setu(s.integrator.sol, s.simple_sys.set_values)
        s.get_pos = getu(s.integrator.sol, s.simple_sys.pos[:,:])
        s.get_flap_angle = getu(s.integrator.sol, s.simple_sys.flap_angle)
        s.get_flap_acc = getu(s.integrator.sol, s.simple_sys.flap_acc)
        s.get_vel_kite = getu(s.integrator.sol, s.simple_sys.vel[:,s.num_A])
        s.get_winch_forces = getu(s.integrator.sol, s.simple_sys.winch_force)
        s.get_L_C = getu(s.integrator.sol, s.simple_sys.L_C)
        s.get_L_D = getu(s.integrator.sol, s.simple_sys.L_D)
        s.get_D_C = getu(s.integrator.sol, s.simple_sys.D_C)
        s.get_D_D = getu(s.integrator.sol, s.simple_sys.D_D)
        s.get_tether_lengths = getu(s.integrator.sol, s.simple_sys.tether_length)
        s.get_tether_vels = getu(s.integrator.sol, s.simple_sys.tether_vel)
        s.get_heading = getu(s.integrator.sol, s.simple_sys.heading)
    end
    s.set_set_values(s.integrator, set_values)
    s.integrator.p[s.v_wind_gnd_idx] .= s.v_wind_gnd
    s.integrator.p[s.v_wind_idx] .= s.v_wind
    s.t_0 = s.integrator.t
    OrdinaryDiffEqCore.step!(s.integrator, dt, true)
    if !successful_retcode(s.integrator.sol)
        println("Return code for solution: ", s.integrator.sol.retcode)
    end
    @assert successful_retcode(s.integrator.sol)
    update_state!(s)
    s.integrator.t
end

function calc_pre_tension(s::KPS4_3L)
    forces = spring_forces(s)
    avg_force = 0.0
    for i in 1:s.num_A
        avg_force += forces[i]
    end
    avg_force /= s.num_A
    res = avg_force/s.c_spring
    if res < 0.0 res = 0.0 end
    if isnan(res) res = 0.0 end
    return res + 1.0
end

"""
    unstretched_length(s::KPS4_3L)

Getter for the unstretched tether reel-out lenght (at zero force).
"""
function unstretched_length(s::KPS4_3L) s.tether_lengths[3] end

"""
    tether_length(s::KPS4_3L)

Calculate and return the real, stretched tether lenght.
"""
function tether_length(s::KPS4_3L)
    length = 0.0
    for i in 3:3:s.num_flap_C-1
        length += norm(s.pos[i+3] - s.pos[i])
    end
    return length
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
    s.pos[end]
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
    winch_force(s::KPS4_3L)

Return the absolute value of the force at the winch as calculated during the last timestep. 
"""
function winch_force(s::KPS4_3L) s.get_winch_forces(s.integrator) end


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
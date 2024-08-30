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


const KITE_SPRINGS_3L = 6
const KITE_PARTICLES_3L = 4
const KITE_ANGLE_3L = 0.0

"""
    mutable struct KPS4_3L{S, T, P, Q, SP} <: AbstractKiteModel

State of the kite power system, using a 3 point kite model and three steering lines to the ground. Parameters:
- S: Scalar type, e.g. SimFloat
  In the documentation mentioned as Any, but when used in this module it is always SimFloat and not Any.
- T: Vector type, e.g. MVector{3, SimFloat}
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
    last_init_elevation::SimFloat     = 0.0
    "The last initial tether length"
    last_init_tether_length::SimFloat = 0.0
    "Reference to the last settings hash"
    last_set_hash::UInt64   = 0
    "Reference to the atmospheric model as implemented in the package AtmosphericModels"
    am::AtmosphericModel = AtmosphericModel()
    "Function for calculation the lift coefficent, using a spline based on the provided value pairs."
    calc_cl = Spline1D(se().alpha_cl, se().cl_list)
    "Function for calculation the drag coefficent, using a spline based on the provided value pairs."
    calc_cd = Spline1D(se().alpha_cd, se().cd_list)
    "Reference to the motor models as implemented in the package WinchModels. index 1: middle motor, index 2: left motor, index 3: right motor"
    motors::SVector{3, AbstractWinchModel}
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
    "last winch force"
    winch_forces::SVector{3,T} = [zeros(S, 3) for _ in 1:3]
    "a copy of the residual one (pos,vel) for debugging and unit tests"    
    res1::SVector{P, T} = zeros(SVector{P, T})
    "a copy of the residual two (vel,acc) for debugging and unit tests"
    res2::SVector{P, T} = zeros(SVector{P, T})
    "a copy of the actual positions as output for the user"
    pos::SVector{P, T} = zeros(SVector{P, T})
    vel::SVector{P, T} = zeros(SVector{P, T})
    veld::SVector{P, T} = zeros(SVector{P, T})
    "velocity vector of the kite"
    vel_kite::T =          zeros(S, 3)
    "unstressed segment lengths of the three tethers [m]"
    segment_lengths::T =           zeros(S, 3)
    "azimuth angle in radian; inital value is zero"
    psi::S =              zero(S)
    "relative start time of the current time interval"
    t_0::S =               0.0
    "reel out speed of the winch"
    reel_out_speeds::T =        zeros(S, 3)
    "unstretched tether length"
    tether_lengths::T =          zeros(S, 3)
    "lengths of the connections of the steering tethers to the kite"
    steering_pos::MVector{2, S} =      zeros(S, 2)
    "air density at the height of the kite"
    rho::S =               0.0
    "multiplier for the stiffniss of tether and bridle"
    stiffness_factor::S =  1.0
    "initial masses of the point masses"
    initial_masses::MVector{P, S} = ones(P)
    "current masses, depending on the total tether length"
    masses::MVector{P, S}         = zeros(P)
    "vector of the springs, defined as struct"
    springs::MVector{Q, SP}       = zeros(SP, Q)
    "vector of the forces, acting on the particles"
    forces::SVector{P, T} = zeros(SVector{P, T})
    "synchronous speed or torque of the motor/ generator"
    set_values::KVec3  = zeros(KVec3)
    torque_control::Bool = false
    "x vector of kite reference frame"
    e_x::T =                 zeros(S, 3)
    "y vector of kite reference frame"
    e_y::T =                 zeros(S, 3)
    "z vector of kite reference frame"
    e_z::T =                 zeros(S, 3)
    e_r::T =                 zeros(S, 3)
    "Point number of E"
    num_E::Int64 =           0
    "Point number of C"
    num_C::Int64 =           0
    "Point number of D"
    num_D::Int64 =           0
    "Point number of A"
    num_A::Int64 =           0
    "Angle of right tip"
    α_l::SimFloat =     0.0
    "Angle of left tip"
    α_r::SimFloat =     0.0

    L_C::T = zeros(S, 3)
    L_D::T = zeros(S, 3)
    D_C::T = zeros(S, 3)
    D_D::T = zeros(S, 3)

    steady_sol::Union{SciMLBase.NonlinearSolution, Nothing} = nothing
    simple_sys::Union{ModelingToolkit.ODESystem, Nothing} = nothing

    set_values_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    v_wind_gnd_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    stiffness_factor_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    v_wind_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    prob::Union{OrdinaryDiffEqCore.ODEProblem, Nothing} = nothing
    get_pos::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_steering_pos::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_line_acc::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_kite_vel::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_winch_forces::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_tether_lengths::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_tether_speeds::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_L_C::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_L_D::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_D_C::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_D_D::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing

    # half_drag_force::SVector{P, T} = zeros(SVector{P, T})
end

"""
    clear!(s::KPS4_3L)

Initialize the kite power model.
"""
function clear!(s::KPS4_3L)
    s.iter = 0
    s.t_0 = 0.0                              # relative start time of the current time interval
    s.reel_out_speeds = zeros(3)
    # s.last_reel_out_speeds = zeros(3)
    s.v_wind_gnd    .= [s.set.v_wind, 0.0, 0.0]    # wind vector at reference height
    s.v_wind_tether .= [s.set.v_wind, 0.0, 0.0]
    s.v_apparent    .= [s.set.v_wind, 0.0, 0.0]
    height = sin(deg2rad(s.set.elevation)) * (s.set.l_tether)
    s.v_wind .= s.v_wind_gnd * calc_wind_factor(s.am, height)
    s.L_C .= 0.0
    s.L_D .= 0.0
    s.D_C .= 0.0
    s.D_D .= 0.0
    s.e_x .= 0.0
    s.e_y .= 0.0
    s.e_z .= 0.0
    s.set_values .= 0.0
    s.steering_pos .= 0.0
    s.vel_kite .= 0.0
    [s.winch_forces[i] .= 0.0 for i in 1:3]
    s.tether_lengths .= [s.set.l_tether for _ in 1:3]
    s.α_l = π/2 - s.set.min_steering_line_distance/(2*s.set.radius)
    s.α_r = π/2 + s.set.min_steering_line_distance/(2*s.set.radius)
    s.segment_lengths .= s.tether_lengths ./ s.set.segments
    s.num_E = s.set.segments*3+3
    s.num_C = s.set.segments*3+3+1
    s.num_D = s.set.segments*3+3+2
    s.num_A = s.set.segments*3+3+3
    for i in 1:s.num_A
        s.forces[i] .= 0.0
        s.veld[i] .= 0.0
    end
    s.rho = s.set.rho_0
    init_masses!(s)
    init_springs!(s)
    s.calc_cl = Spline1D(s.set.alpha_cl, s.set.cl_list)
    s.calc_cd = Spline1D(s.set.alpha_cd, s.set.cd_list) 
end

function KPS4_3L(kcu::KCU)
    set = kcu.set
    if set.winch_model == "TorqueControlledMachine"
        s = KPS4_3L{SimFloat, KVec3, set.segments*3+2+KITE_PARTICLES, set.segments*3+KITE_SPRINGS_3L, SP}(set=kcu.set, motors=[TorqueControlledMachine(set) for _ in 1:3])
    else
        s = KPS4_3L{SimFloat, KVec3, set.segments*3+2+KITE_PARTICLES, set.segments*3+KITE_SPRINGS_3L, SP}(set=kcu.set, motors=[AsyncMachine(set) for _ in 1:3])
    end
    s.num_E = s.set.segments*3+3
    s.num_C = s.set.segments*3+3+1
    s.num_D = s.set.segments*3+3+2
    s.num_A = s.set.segments*3+3+3     
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
    ss.heading = calc_heading(s)
    ss.course = calc_course(s)
    ss.v_app = norm(s.v_apparent)
    ss.l_tether = s.tether_lengths[3]
    ss.v_reelout = s.reel_out_speeds[3]
    ss.depower = 100 - ((s.steering_pos[1] + s.steering_pos[2])/2) / ((s.set.middle_length + s.set.tip_length)/2) * 100
    ss.steering = (s.steering_pos[2] - s.steering_pos[1]) / ((s.set.middle_length + s.set.tip_length)/2) * 100
    ss.vel_kite .= s.vel_kite
    nothing
end

function SysState(s::KPS4_3L, zoom=1.0)
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
    heading = calc_heading(s)
    course = calc_course(s)
    v_app_norm = norm(s.v_apparent)
    t_sim = 0
    depower = 100 - ((s.steering_pos[1] + s.steering_pos[2])/2) / ((s.set.middle_length + s.set.tip_length)/2) * 100
    steering = (s.steering_pos[1] - s.steering_pos[2]) / ((s.set.middle_length + s.set.tip_length)/2) * 100
    KiteUtils.SysState{P}(s.t_0, t_sim, 0, 0, orient, elevation, azimuth, s.tether_lengths[3], s.reel_out_speeds[3], forces[3], depower, steering, 
                          heading, course, v_app_norm, s.vel_kite, X, Y, Z, 
                          0, 0, 0, 0, 
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end


"""
    init_sim!(s; t_end=1.0, stiffness_factor=1.0, prn=false)

Initialises the integrator of the model.

Parameters:
- s:     an instance of an abstract kite model
- t_end: end time of the simulation; normally not needed
- stiffness_factor: factor applied to the tether stiffness during initialisation
- prn: if set to true, print the detailed solver results
- steady_state_history: an instance of SteadyStateHistory containing old pairs of AKM objects and integrators

Returns:
An instance of a DAE integrator.
"""
function init_sim!(s::KPS4_3L; t_end=1.0, stiffness_factor=1.0, prn=false, 
                   torque_control=true)
    clear!(s)
    change_control_mode = s.torque_control != torque_control
    s.stiffness_factor = stiffness_factor
    s.torque_control = torque_control
    dt = 1/s.set.sample_freq
    tspan   = (0.0, dt) 
    solver = KenCarp4(autodiff=false) # TRBDF2, Rodas4P, Rodas5P, Kvaerno5, KenCarp4, radau, QNDF

    new_inital_conditions = (s.last_init_elevation != s.set.elevation || s.last_init_tether_length != s.set.l_tether)
    s.set_hash = settings_hash(s.set)
    steady_tol = 1e-6
    if isnothing(s.prob) || change_control_mode || s.last_set_hash != s.set_hash
        if prn; println("initializing with new model and new steady state"); end
        pos = init_pos(s)
        model!(s, pos; torque_control=s.torque_control)
        s.prob = ODEProblem(s.simple_sys, nothing, tspan)
        steady_prob = SteadyStateProblem(s.prob)
        s.steady_sol = solve(steady_prob, DynamicSS(solver; tspan=tspan), abstol=steady_tol, reltol=steady_tol)
        s.prob = remake(s.prob; u0=s.steady_sol.u)
    elseif new_inital_conditions
        if prn; println("initializing with last model and new steady state"); end
        pos = init_pos(s)
        s.prob = ODEProblem(s.simple_sys, [s.simple_sys.pos => pos, s.simple_sys.tether_length => s.tether_lengths], tspan)
        steady_prob = SteadyStateProblem(s.prob)
        s.steady_sol = solve(steady_prob, DynamicSS(solver; tspan=tspan), abstol=steady_tol, reltol=steady_tol)
        s.prob = remake(s.prob; u0=s.steady_sol.u)
    else
        if prn; println("initializing with last model and last steady state"); end
    end
    s.last_init_elevation = deepcopy(s.set.elevation)
    s.last_init_tether_length = deepcopy(s.set.l_tether)    
    s.last_set_hash = deepcopy(s.set_hash)
    integrator = OrdinaryDiffEqCore.init(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
    if isnothing(s.set_values_idx)
        s.set_values_idx = parameter_index(integrator.f, :set_values)
        s.v_wind_gnd_idx = parameter_index(integrator.f, :v_wind_gnd)
        s.v_wind_idx = parameter_index(integrator.f, :v_wind)
        s.stiffness_factor_idx = parameter_index(integrator.f, :stiffness_factor)
        s.get_pos = getu(integrator.sol, s.simple_sys.pos[:,:])
        s.get_steering_pos = getu(integrator.sol, s.simple_sys.steering_pos)
        s.get_line_acc = getu(integrator.sol, s.simple_sys.acc[:,s.num_E-2])
        s.get_kite_vel = getu(integrator.sol, s.simple_sys.vel[:,s.num_A])
        s.get_winch_forces = getu(integrator.sol, s.simple_sys.force[:,1:3])
        s.get_L_C = getu(integrator.sol, s.simple_sys.L_C)
        s.get_L_D = getu(integrator.sol, s.simple_sys.L_D)
        s.get_D_C = getu(integrator.sol, s.simple_sys.D_C)
        s.get_D_D = getu(integrator.sol, s.simple_sys.D_D)
        s.get_tether_lengths = getu(integrator.sol, s.simple_sys.tether_length)
        s.get_tether_speeds = getu(integrator.sol, s.simple_sys.tether_speed)
    end
    update_pos!(s, integrator)
    return integrator
end

function next_step!(s::KPS4_3L, integrator; set_values=zeros(KVec3), v_wind_gnd=s.set.v_wind, wind_dir=0.0, dt=1/s.set.sample_freq)
    s.iter = 0
    set_v_wind_ground!(s, calc_height(s), v_wind_gnd, wind_dir)
    s.set_values .= set_values
    integrator.ps[s.set_values_idx] .= s.set_values
    integrator.ps[s.v_wind_gnd_idx] .= s.v_wind_gnd
    integrator.ps[s.v_wind_idx] .= s.v_wind
    integrator.ps[s.stiffness_factor_idx] = s.stiffness_factor
    s.t_0 = integrator.t
    OrdinaryDiffEqCore.step!(integrator, dt, true)
    update_pos!(s, integrator)
    if s.stiffness_factor < 1.0
        s.stiffness_factor+=0.01
        if s.stiffness_factor > 1.0
            s.stiffness_factor = 1.0
        end
    end
    integrator.t
end

function calc_pre_tension(s::KPS4_3L)
    forces = spring_forces(s)
    avg_force = 0.0
    for i in 1:s.num_A
        avg_force += forces[i]
    end
    avg_force /= s.num_A
    res = avg_force/s.set.c_spring
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
    for i in 3:3:s.num_E-3
        length += norm(s.pos[i+3] - s.pos[i])
    end
    return length
end



# =================== getter functions ====================================================

"""
    calc_height(s::KPS4_3L)

Determine the height of the topmost kite particle above ground.
"""
function calc_height(s::KPS4_3L)isnothing(s.prob) || 
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
    kite_ref_frame(s::KPS4_3L)

Returns a tuple of the x, y, and z vectors of the kite reference frame.
"""
function kite_ref_frame(s::KPS4_3L)
    s.e_x, s.e_y, s.e_z
end

"""
    winch_force(s::KPS4_3L)

Return the absolute value of the force at the winch as calculated during the last timestep. 
"""
function winch_force(s::KPS4_3L) norm.(s.winch_forces) end


# ==================== mtk model functions ================================================
# Implementation of the three-line model using ModellingToolkit.jl


# issue: still uses settings getter function
function calc_acc_speed(tether_speed::SimFloat, norm_::SimFloat, set_speed::SimFloat)
    calc_acceleration(AsyncMachine(se("system_3l.yaml")), tether_speed, norm_; set_speed, set_torque=nothing, use_brake=false)
end
@register_symbolic calc_acc_speed(tether_speed, norm_, set_speed)

function calc_acc_torque(tether_speed::SimFloat, norm_::SimFloat, set_torque::SimFloat)
    calc_acceleration(TorqueControlledMachine(se("system_3l.yaml")), tether_speed, norm_; set_speed=nothing, set_torque, use_brake=false)
end
@register_symbolic calc_acc_torque(tether_speed, norm_, set_torque)

"""
    calc_aero_forces!(s::KPS4_3L, eqs2, force_eqs, force, pos, vel, t, e_x, e_y, e_z, rho)

Calculates the aerodynamic forces acting on the kite particles.

Parameters:
- pos:              vector of the particle positions
- vel:              vector of the particle velocities
- rho:              air density [kg/m^3]

Updates the vector s.forces of the first parameter.
"""
function calc_aero_forces_mtk!(s::KPS4_3L, eqs2, force_eqs, force, pos, vel, t, e_x, e_y, e_z, rho, v_wind, steering_pos)
    n = s.set.aero_surfaces
    @variables begin
        E_c(t)[1:3]
        v_cx(t)[1:3]
        v_dx(t)[1:3]
        v_dy(t)[1:3]
        v_dz(t)[1:3]
        v_cy(t)[1:3]
        v_cz(t)[1:3]
        y_lc(t)
        y_ld(t)
    end

    eqs2 = [
        eqs2
        # in the aero calculations, E_c is the center of the circle shape on which the kite lies
        E_c     ~ pos[:, s.num_E] + e_z * (-s.set.bridle_center_distance + s.set.radius)
        v_cx    ~ (vel[:, s.num_C] ⋅ e_x) * e_x
        v_dx    ~ (vel[:, s.num_D] ⋅ e_x) * e_x
        v_dy    ~ (vel[:, s.num_D] ⋅ e_y) * e_y
        v_dz    ~ (vel[:, s.num_D] ⋅ e_z) * e_z
        v_cy    ~ (vel[:, s.num_C] ⋅ e_y) * e_y
        v_cz    ~ (vel[:, s.num_C] ⋅ e_z) * e_z
        y_lc    ~  norm(pos[:, s.num_C] - 0.5 * (pos[:, s.num_C] + pos[:, s.num_D]))
        y_ld    ~ -norm(pos[:, s.num_D] - 0.5 * (pos[:, s.num_C] + pos[:, s.num_D]))
    ]

    # integrating loop variables
    @variables begin
        F(t)[1:3, 1:2n]
        e_r(t)[1:3, 1:2n]
        y_l(t)[1:2n]
        v_kite(t)[1:3, 1:2n]
        v_a(t)[1:3, 1:2n]
        e_drift(t)[1:3, 1:2n]
        v_a_xr(t)[1:3, 1:2n]
        aoa(t)[1:n*2]
        dL_dα(t)[1:3, 1:2n]
        dD_dα(t)[1:3, 1:2n]
        L_C(t)[1:3]
        L_D(t)[1:3]
        D_C(t)[1:3]
        D_D(t)[1:3]
        F_steering_c(t)[1:3]
        F_steering_d(t)[1:3]
        d(t)[1:2n]
    end
    l_c_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(L_C .~ 0))
    l_d_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(L_D .~ 0))
    d_c_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(D_C .~ 0))
    d_d_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(D_D .~ 0))
    kite_length = zeros(MVector{2n, SimFloat})
    α           = zero(SimFloat)
    α_0         = zero(SimFloat)
    α_middle    = zero(SimFloat)
    dα          = zero(SimFloat)
    # Calculate the lift and drag
    α_0         = π/2 - s.set.width/2/s.set.radius
    α_middle    = π/2
    dα          = (α_middle - α_0) / n
    for i in 1:n*2
        if i <= n
            α = α_0 + -dα/2 + i * dα
        else
            α = pi - (α_0 + -dα/2 + (i-n) * dα)
        end

        eqs2 = [
            eqs2
            F[:, i]          ~ E_c + e_y * cos(α) * s.set.radius - e_z * sin(α) * s.set.radius
            e_r[:, i]        ~ (E_c - F[:, i]) / norm(E_c - F[:, i])
            y_l[i]           ~ cos(α) * s.set.radius
            α < π/2 ?
                v_kite[:, i] ~ ((v_cx - v_dx) / (y_lc - y_ld) * (y_l[i] - y_ld) + v_dx) + v_cy + v_cz :
                v_kite[:, i] ~ ((v_cx - v_dx) / (y_lc - y_ld) * (y_l[i] - y_ld) + v_dx) + v_dy + v_dz
            v_a[:,i]         ~ v_wind .- v_kite[:,i]
            e_drift[:, i]    ~ (e_r[:, i] × e_x)
            v_a_xr[:, i]     ~ v_a[:, i] .- (v_a[:, i] ⋅ e_drift[:, i]) .* e_drift[:, i]
        ]
        if α < π/2
            kite_length[i] = (s.set.tip_length + (s.set.middle_length-s.set.tip_length) * α 
                              * s.set.radius/(0.5*s.set.width))
        else
            kite_length[i] = (s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (π-α) 
                              * s.set.radius/(0.5*s.set.width))
        end
        eqs2 = [
            eqs2
            α < s.α_l ?
                d[i]    ~ steering_pos[1] :
            α > s.α_r ?
                d[i]    ~ steering_pos[2] :
                d[i]    ~ (steering_pos[2] - steering_pos[1]) / (s.α_r - s.α_l) * (α - s.α_l) + (steering_pos[1])
            aoa[i]      ~ -asin((v_a_xr[:, i] / norm(v_a_xr[:, i])) ⋅ e_r[:, i]) + 
                           asin(clamp(d[i] / kite_length[i], -1.0, 1.0))
            dL_dα[:, i] ~ 0.5 * rho * (norm(v_a_xr[:, i]))^2 * s.set.radius * kite_length[i] * rad_cl_mtk(aoa[i]) * 
                                ((v_a_xr[:, i] × e_drift[:, i]) / norm(v_a_xr[:, i] × e_drift[:, i]))
            dD_dα[:, i] ~ 0.5 * rho * norm(v_a_xr[:, i]) * s.set.radius * kite_length[i] * rad_cd_mtk(aoa[i]) * 
                                v_a_xr[:,i] # the sideways drag cannot be calculated with the C_d formula
        ]
        if i <= n
            [l_c_eq[j] = (L_C[j] ~ l_c_eq[j].rhs + dL_dα[j, i] * dα) for j in 1:3]
            [d_c_eq[j] = (D_C[j] ~ d_c_eq[j].rhs + dD_dα[j, i] * dα) for j in 1:3]
        else 
            [l_d_eq[j] = (L_D[j] ~ l_d_eq[j].rhs + dL_dα[j, i] * dα) for j in 1:3]
            [d_d_eq[j] = (D_D[j] ~ d_d_eq[j].rhs + dD_dα[j, i] * dα) for j in 1:3]
        end
    end
    
    eqs2 = [
        eqs2
        l_c_eq
        d_c_eq
        l_d_eq
        d_d_eq
        F_steering_c ~ ((0.2 * (L_C ⋅ -e_z)) * -e_z)
        F_steering_d ~ ((0.2 * (L_D ⋅ -e_z)) * -e_z)
    ]
    
    force_eqs[:,s.num_C]   .= (force[:,s.num_C]   .~ (L_C + D_C) - F_steering_c) 
    force_eqs[:,s.num_D]   .= (force[:,s.num_D]   .~ (L_D + D_D) - F_steering_d)
    force_eqs[:,s.num_E-2] .= (force[:,s.num_E-2] .~ F_steering_c)
    force_eqs[:,s.num_E-1] .= (force[:,s.num_E-1] .~ F_steering_d)
    return eqs2, force_eqs
end

""" 
    calc_particle_forces!(s::KPS4_3L, eqs2, force_eqs, force, pos1, pos2, vel1, vel2, length, c_spring, damping, 
                          rho, i, l_0, k, c, segment, rel_vel, av_vel, norm1, unit_vector, k1, k2, c1, spring_vel,
                          spring_force, v_apparent, v_wind_tether, area, v_app_perp, half_drag_force)

Calculate the drag force and spring force of the tether segment, defined by the parameters pos1, pos2, vel1 and vel2
and distribute it equally on the two particles, that are attached to the segment.
The result is stored in the array s.forces. 
"""
function calc_particle_forces_mtk!(s::KPS4_3L, eqs2, force_eqs, force, pos1, pos2, vel1, vel2, length, c_spring, 
    damping, rho, i, l_0, k, c, segment, rel_vel, av_vel, norm1, unit_vector, k1, k2, c1, spring_vel,
            spring_force, v_apparent, v_wind_tether, area, v_app_perp, half_drag_force, stiffness_factor)
    d_tether = s.set.d_tether/1000.0
    eqs2 = [
        eqs2
        i <= s.set.segments*3 ? l_0 ~ length[(i-1) % 3 + 1] : l_0 ~ s.springs[i].length # Unstressed length
        i <= s.set.segments*3 ? k   ~ c_spring[(i-1) % 3 + 1] * stiffness_factor :
                                k   ~ s.springs[i].c_spring * stiffness_factor        # Spring constant
        i <= s.set.segments*3 ? c   ~ damping[(i-1) % 3 + 1] : c ~ s.springs[i].damping # Damping coefficient    
        segment     .~ pos1 - pos2
        rel_vel     .~ vel1 - vel2
        av_vel      .~ 0.5 * (vel1 + vel2)
        norm1        ~ norm(segment)
        unit_vector .~ segment / norm1
        k1           ~ 0.25 * k # compression stiffness kite segments
        k2           ~ 0.1 * k  # compression stiffness tether segments
        c1           ~ 6.0 * c  # damping kite segments
        spring_vel  .~ rel_vel ⋅ unit_vector
    ]

    if i >= s.num_E-2  # kite springs
        for j in 1:3
            eqs2 = [
                eqs2
                spring_force[j] ~ ifelse(
                    (norm1 - l_0) > 0.0,
                    (k  * (l_0 - norm1) - c1 * spring_vel) * unit_vector[j],
                    (k1 * (l_0 - norm1) -  c * spring_vel) * unit_vector[j]
                )
            ]
        end
    else
        for j in 1:3
            eqs2 = [
                eqs2
                spring_force[j] ~ ifelse(
                    (norm1 - l_0) > 0.0,
                    (k  * (l_0 - norm1) - c * spring_vel) * unit_vector[j],
                    (k2 * (l_0 - norm1) - c * spring_vel) * unit_vector[j]
                    )
            ]
        end
    end
    eqs2 = [
        eqs2
        v_apparent       ~ v_wind_tether - av_vel
        area             ~ norm1 * d_tether
        v_app_perp       ~ v_apparent - v_apparent ⋅ unit_vector * unit_vector
        half_drag_force .~ (0.25 * rho * s.set.cd_tether * norm(v_app_perp) * area) .* v_app_perp
    ]

    for j in 1:3
        force_eqs[j, s.springs[i].p1] = 
            (force[j,s.springs[i].p1] ~ force_eqs[j, s.springs[i].p1].rhs + (half_drag_force[j] + spring_force[j]))
        force_eqs[j, s.springs[i].p2] = 
            (force[j,s.springs[i].p2] ~ force_eqs[j, s.springs[i].p2].rhs + (half_drag_force[j] - spring_force[j]))
    end
    
    return eqs2, force_eqs
end



"""
    inner_loop_mtk!(s::KPS4_3L, eqs2, force_eqs, t, force, pos, vel, length, c_spring, damping, v_wind_gnd)

Calculate the forces, acting on all particles.

Output:length
- s.forces
- s.v_wind_tether
"""
@inline function inner_loop_mtk!(s::KPS4_3L, eqs2, force_eqs, t, force, pos, vel, length, c_spring, damping, v_wind_gnd, stiffness_factor)
    @variables begin
        height(t)[eachindex(s.springs)]
        rho(t)[eachindex(s.springs)]
        v_wind_tether(t)[1:3, eachindex(s.springs)]

        l_0(t)[eachindex(s.springs)]
        k(t)[eachindex(s.springs)]
        c(t)[eachindex(s.springs)]
        segment(t)[1:3,eachindex(s.springs)]
        rel_vel(t)[1:3, eachindex(s.springs)]
        av_vel(t)[1:3, eachindex(s.springs)] 
        norm1(t)[eachindex(s.springs)]
        unit_vector(t)[1:3, eachindex(s.springs)]
        k1(t)[eachindex(s.springs)]
        k2(t)[eachindex(s.springs)]
        c1(t)[eachindex(s.springs)]
        spring_vel(t)[eachindex(s.springs)]
        spring_force(t)[1:3, eachindex(s.springs)]
        v_apparent(t)[1:3, eachindex(s.springs)]
        area(t)[eachindex(s.springs)]
        v_app_perp(t)[1:3, eachindex(s.springs)]
        half_drag_force(t)[1:3, eachindex(s.springs)]
    end
    
    for i in eachindex(s.springs)
        p1 = s.springs[i].p1  # First point nr.
        p2 = s.springs[i].p2
        eqs2 = [
            eqs2
            height[i]           ~ 0.5 * (pos[:, p1][3] + pos[:, p2][3])
            rho[i]              ~ calc_rho(s.am, height[i])
            v_wind_tether[:, i] ~ calc_wind_factor(s.am, height[i]) * v_wind_gnd
        ]

        # TODO: @assert height > 0
        eqs2, force_eqs = calc_particle_forces_mtk!(s, eqs2, force_eqs, force, pos[:, p1], pos[:, p2], vel[:, p1], 
                          vel[:, p2], length, c_spring, damping, rho[i], i, l_0[i], k[i], c[i], segment[:, i], 
                          rel_vel[:, i], av_vel[:, i], norm1[i], unit_vector[:, i], k1[i], k2[i], c1[i], spring_vel[i],
                          spring_force[:, i], v_apparent[:,i], v_wind_tether[:, i], area[i], v_app_perp[:, i], 
                          half_drag_force[:, i], stiffness_factor)
    end

    return eqs2, force_eqs
end

function update_pos!(s, integrator)
    pos = s.get_pos(integrator)
    s.steering_pos     .= s.get_steering_pos(integrator)
    [s.pos[i]          .= pos[:,i] for i in 1:s.num_A]
    s.veld[s.num_E-2]  .= s.get_line_acc(integrator)
    s.vel_kite         .= s.get_kite_vel(integrator)
    winch_forces        = s.get_winch_forces(integrator)
    [s.winch_forces[i] .= (winch_forces[:,i]) for i in 1:3]
    s.tether_lengths   .= s.get_tether_lengths(integrator)
    s.reel_out_speeds  .= s.get_tether_speeds(integrator)
    s.L_C               = s.get_L_C(integrator)
    s.L_D               = s.get_L_D(integrator)
    s.D_C               = s.get_D_C(integrator)
    s.D_D               = s.get_D_D(integrator)
    calc_kite_ref_frame!(s, s.pos[s.num_E], s.pos[s.num_C], s.pos[s.num_D])
    @assert all(abs.(s.steering_pos) .<= s.set.tip_length)
    nothing
end

function model!(s::KPS4_3L, pos_; torque_control=false)
    @parameters begin
        set_values[1:3] = s.set_values
        v_wind_gnd[1:3] = s.v_wind_gnd
        v_wind[1:3] = s.v_wind
        stiffness_factor = s.stiffness_factor
    end
    @variables begin
        pos(t)[1:3, 1:s.num_A] = pos_ # left right middle
        vel(t)[1:3, 1:s.num_A] = zeros(3, s.num_A) # left right middle
        acc(t)[1:3, 1:s.num_A] = zeros(3, s.num_A) # left right middle
        tether_length(t)[1:3]  = s.tether_lengths
        steering_pos(t)[1:2]   = s.steering_pos
        steering_vel(t)[1:2]   = zeros(2)
        steering_acc(t)[1:2]   = zeros(2)
        tether_speed(t)[1:3]   = zeros(3) # left right middle
        segment_length(t)[1:3] = zeros(3) # left right middle
        mass_tether_particle(t)[1:3]      # left right middle
        damping(t)[1:3] = s.set.damping ./ s.tether_lengths ./ s.set.segments   # left right middle
        c_spring(t)[1:3] = s.set.c_spring ./ s.tether_lengths ./ s.set.segments # left right middle
        P_c(t)[1:3] = 0.5 .* (s.pos[s.num_C] + s.pos[s.num_D])
        e_x(t)[1:3]
        e_y(t)[1:3]
        e_z(t)[1:3]
        force(t)[1:3, 1:s.num_A] = zeros(3, s.num_A) # left right middle
        rho_kite(t) = 0.0
    end
    # Collect the arrays into variables
    pos = collect(pos)
    vel = collect(vel)
    acc = collect(acc)

    eqs1 = []
    mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2

    [eqs1 = vcat(eqs1, D.(pos[:, i]) .~ 0.0) for i in 1:3]
    [eqs1 = vcat(eqs1, D.(pos[:, i]) .~ vel[:,i]) for i in 4:s.num_E-3]
    eqs1 = [eqs1; D.(steering_pos)   .~ steering_vel]
    [eqs1 = vcat(eqs1, D.(pos[:, i]) .~ vel[:,i]) for i in s.num_E:s.num_A]
    [eqs1 = vcat(eqs1, D.(vel[:, i]) .~ 0.0) for i in 1:3]
    [eqs1 = vcat(eqs1, D.(vel[:, i]) .~ acc[:,i]) for i in 4:s.num_E-3]
    eqs1 = [eqs1; D.(steering_vel)   .~ steering_acc]
    [eqs1 = vcat(eqs1, D.(vel[:, i]) .~ acc[:,i]) for i in s.num_E:s.num_A]

    eqs1 = vcat(eqs1, D.(tether_length) .~ tether_speed)
    if torque_control
        eqs1 = vcat(eqs1, D.(tether_speed) .~ [calc_acc_torque(tether_speed[i], norm(force[:, (i-1) % 3 + 1]),
                                                               set_values[i]) for i in 1:3])
    else
        eqs1 = vcat(eqs1, D.(tether_speed) .~ [calc_acc_speed(tether_speed[i], norm(force[:,(i-1) % 3 + 1]), 
                                                               set_values[i]) for i in 1:3])
    end

    # Compute the masses and forces
    force_eqs = SizedArray{Tuple{3, s.num_A}, Symbolics.Equation}(undef)
    force_eqs[:, :] .= (force[:, :] .~ 0)

    eqs2 = [
        pos[:, s.num_E-2] ~ pos[:, s.num_C] + e_z * steering_pos[1]
        pos[:, s.num_E-1] ~ pos[:, s.num_D] + e_z * steering_pos[2]
        vel[:, s.num_E-2] ~ vel[:, s.num_C] + e_z * steering_vel[1]
        vel[:, s.num_E-1] ~ vel[:, s.num_D] + e_z * steering_vel[2]
        acc[:, s.num_E-2] ~ acc[:, s.num_C] + e_z * steering_acc[1]
        acc[:, s.num_E-1] ~ acc[:, s.num_D] + e_z * steering_acc[2]
        segment_length       ~ tether_length  ./ s.set.segments
        mass_tether_particle ~ mass_per_meter .* segment_length
        damping              ~ s.set.damping  ./ segment_length
        c_spring             ~ s.set.c_spring ./ segment_length
        P_c ~ 0.5 * (pos[:, s.num_C] + pos[:, s.num_D])
        e_y ~ (pos[:, s.num_C] - pos[:, s.num_D]) / norm(pos[:, s.num_C] - pos[:, s.num_D])
        e_z ~ (pos[:, s.num_E] - P_c) / norm(pos[:, s.num_E] - P_c)
        e_x ~ cross(e_y, e_z)
        rho_kite ~ calc_rho(s.am, pos[3,s.num_A])
    ]

    eqs2, force_eqs = calc_aero_forces_mtk!(s, eqs2, force_eqs, force, pos, vel, t, e_x, e_y, e_z, rho_kite, v_wind, steering_pos)
    eqs2, force_eqs = inner_loop_mtk!(s, eqs2, force_eqs, t, force, pos, vel, segment_length, c_spring, damping, 
                                      v_wind_gnd, stiffness_factor)
    
    for i in 1:3
        eqs2 = vcat(eqs2, vcat(force_eqs[:, i]))
        eqs2 = vcat(eqs2, acc[:, i] .~ 0)
    end
    for i in 4:s.num_E-3
        eqs2 = vcat(eqs2, vcat(force_eqs[:,i]))
        eqs2 = vcat(eqs2, acc[:, i] .~ [0.0; 0.0; -G_EARTH] .+ (force[:,i] ./ mass_tether_particle[(i-1)%3+1]))
    end
    for i in s.num_E-2:s.num_E-1
        flap_resistance   = [50.0 * ((vel[:,i]-vel[:, s.num_C]) ⋅ e_z) * e_z[j] for j in 1:3]
        [force_eqs[j,i]   = force[j,i] ~ force_eqs[j,i].rhs + [0.0; 0.0; -G_EARTH][j] + flap_resistance[j] for j in 1:3]
        tether_rhs        = [force_eqs[j, i].rhs for j in 1:3]
        kite_rhs          = [force_eqs[j, i+3].rhs for j in 1:3]
        f_xy              = (tether_rhs ⋅ e_z) .* e_z
        force_eqs[:,i]   .= force[:, i] .~ tether_rhs .- f_xy
        force_eqs[:,i+3] .= force[:, i+3] .~ kite_rhs .+ f_xy
        eqs2              = vcat(eqs2, vcat(force_eqs[:, i]))
        eqs2              = vcat(eqs2, steering_acc[i-s.num_E+3] ~ (force[:,i] ./ mass_tether_particle[(i-1) % 3 + 1]) ⋅ 
                                                                    e_z - (acc[:, i+3] ⋅ e_z))
    end
    for i in s.num_E:s.num_A
        eqs2 = vcat(eqs2, vcat(force_eqs[:,i]))
        eqs2 = vcat(eqs2, acc[:, i] .~ [0.0; 0.0; -G_EARTH] .+ (force[:, i] ./ s.masses[i]))
    end

    eqs = vcat(eqs1, eqs2)

    @named sys = ODESystem(Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs))), t)
    s.simple_sys = structural_simplify(sys)
    nothing
end


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
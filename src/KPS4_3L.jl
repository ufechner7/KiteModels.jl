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
    "drag force of kite and bridle; output of calc_aero_forces!"
    drag_force::T =       zeros(S, 3)
    "lift force of the kite; output of calc_aero_forces!"
    lift_force::T =       zeros(S, 3)    
    "spring force of the current tether segment, output of calc_particle_forces!"
    spring_force::T =     zeros(S, 3)
    "last winch force"
    winch_forces::SVector{3,T} = [zeros(S, 3) for _ in 1:3]
    "a copy of the residual one (pos,vel) for debugging and unit tests"    
    res1::SVector{P, T} = zeros(SVector{P, T})
    "a copy of the residual two (vel,acc) for debugging and unit tests"
    res2::SVector{P, T} = zeros(SVector{P, T})
    "a copy of the actual positions as output for the user"
    pos::SVector{P, T} = zeros(SVector{P, T})
    stable_pos::SVector{P, T} = zeros(SVector{P, T})
    vel::SVector{P, T} = zeros(SVector{P, T})
    posd::SVector{P, T} = zeros(SVector{P, T})
    veld::SVector{P, T} = zeros(SVector{P, T})
    "velocity vector of the kite"
    vel_kite::T =          zeros(S, 3)
    steering_vel::T =          zeros(S, 3)
    "unstressed segment lengths of the three tethers [m]"
    segment_lengths::T =           zeros(S, 3)
    "lift coefficient of the kite, depending on the angle of attack"
    param_cl::S =         0.2
    "drag coefficient of the kite, depending on the angle of attack"
    param_cd::S =         1.0
    "azimuth angle in radian; inital value is zero"
    psi::S =              zero(S)
    "relative start time of the current time interval"
    t_0::S =               0.0
    "reel out speed of the winch"
    reel_out_speeds::T =        zeros(S, 3)
    # "reel out speed at the last time step"
    # last_reel_out_speeds::T =   zeros(S, 3)
    "unstretched tether length"
    tether_lengths::T =          zeros(S, 3)
    "lengths of the connections of the steering tethers to the kite"
    steering_pos::MVector{2, S} =      zeros(S, 2)
    "air density at the height of the kite"
    rho::S =               0.0
    # "actual relative depower setting,  must be between    0 .. 1.0"
    # depower::S =           0.0
    # "actual relative steering setting, must be between -1.0 .. 1.0"
    # steering::S =          0.0
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

    "mtk variables"
    mtk::Bool = false

    set_values_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    v_wind_gnd_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    stiffness_factor_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    v_wind_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    prob::Union{OrdinaryDiffEq.ODEProblem, Nothing} = nothing
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

    half_drag_force::SVector{P, T} = zeros(SVector{P, T})

    "residual variables"
    num_A::Int64 =           0
    L_C::T = zeros(S, 3)
    L_D::T = zeros(S, 3)
    D_C::T = zeros(S, 3)
    D_D::T = zeros(S, 3)
    F_steering_c::T = zeros(S, 3)
    F_steering_d::T = zeros(S, 3)
    dL_dα::T = zeros(S, 3)
    dD_dα::T = zeros(S, 3)
    v_cx::T = zeros(S, 3)
    v_dx::T = zeros(S, 3)
    v_dy::T = zeros(S, 3)
    v_dz::T = zeros(S, 3)
    v_cy::T = zeros(S, 3)
    v_cz::T = zeros(S, 3)
    v_kite::T = zeros(S, 3)
    v_a::T = zeros(S, 3)
    e_drift::T = zeros(S, 3)
    v_a_xr::T = zeros(S, 3)
    E_c::T = zeros(S, 3)
    F::T = zeros(S, 3)
    y_lc::S = 0.0
    y_ld::S = 0.0
    δ_left::S = 0.0
    δ_right::S = 0.0
    α_l::S = 0.0
    α_r::S = 0.0
    distance_c_e::S = 0
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
    s.v_wind_gnd    .= [s.set.v_wind, 0, 0]    # wind vector at reference height
    s.v_wind_tether .= [s.set.v_wind, 0, 0]
    s.v_apparent    .= [s.set.v_wind, 0, 0]
    height = sin(deg2rad(s.set.elevation)) * (s.set.l_tether)
    s.v_wind .= s.v_wind_gnd * calc_wind_factor(s.am, height)

    s.tether_lengths .= [s.set.l_tether for _ in 1:3]
    s.α_l = π/2 - s.set.min_steering_line_distance/(2*s.set.radius)
    s.α_r = π/2 + s.set.min_steering_line_distance/(2*s.set.radius)

    s.segment_lengths .= s.tether_lengths ./ s.set.segments
    s.num_E = s.set.segments*3+3
    s.num_C = s.set.segments*3+3+1
    s.num_D = s.set.segments*3+3+2
    s.num_A = s.set.segments*3+3+3
    init_masses!(s)
    init_springs!(s)
    for i in 1:s.num_A
        s.forces[i] .= zeros(3)
    end
    s.drag_force .= [0.0, 0, 0]
    s.lift_force .= [0.0, 0, 0]
    s.rho = s.set.rho_0
    s.calc_cl = Spline1D(s.set.alpha_cl, s.set.cl_list)
    s.calc_cd = Spline1D(s.set.alpha_cd, s.set.cd_list) 
end


function KPS4_3L(kcu::KCU)
    set = kcu.set
    if set.winch_model == "TorqueControlledMachine"
        s = KPS4_3L{SimFloat, KVec3, set.segments*3+2+KITE_PARTICLES, set.segments*3+KITE_SPRINGS_3L, SP}(set=kcu.set, motors=[TorqueControlledMachine(set) for _ in 1:3])
        println("Using torque control.")
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
    ss.depower = 100 - ((s.δ_left + s.δ_right)/2) / ((s.set.middle_length + s.set.tip_length)/2) * 100
    ss.steering = (s.δ_right - s.δ_left) / ((s.set.middle_length + s.set.tip_length)/2) * 100
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
    depower = 100 - ((s.δ_left + s.δ_right)/2) / ((s.set.middle_length + s.set.tip_length)/2) * 100
    steering = (s.δ_right - s.δ_left) / ((s.set.middle_length + s.set.tip_length)/2) * 100
    KiteUtils.SysState{P}(s.t_0, t_sim, 0, 0, orient, elevation, azimuth, s.tether_lengths[3], s.reel_out_speeds[3], forces[3], depower, steering, 
                          heading, course, v_app_norm, s.vel_kite, X, Y, Z, 
                          0, 0, 0, 0, 
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end


"""
    init_sim!(s; t_end=1.0, stiffness_factor=0.035, prn=false)

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
    s.stiffness_factor = stiffness_factor
    s.torque_control = torque_control
    dt = 1/s.set.sample_freq
    tspan   = (0.0, dt) 
    
    pos = init_pos(s)
    simple_sys, _ = model!(s, pos; torque_control=s.torque_control)
    println("making steady state prob")
    s.prob = ODEProblem(simple_sys, nothing, tspan)
    @time steady_prob = SteadyStateProblem(s.prob)
    println("solving steady state prob")
    @time steady_sol = solve(steady_prob, DynamicSS(KenCarp4(autodiff=false); tspan=tspan), abstol=s.set.abs_tol, reltol=s.set.rel_tol)
    
    solver = KenCarp4(autodiff=false) # TRBDF2, Rodas4P, Rodas5P, Kvaerno5, KenCarp4, radau, QNDF

    s.prob = remake(s.prob; u0=steady_sol.u)
    integrator = OrdinaryDiffEq.init(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
    s.set_values_idx = parameter_index(integrator.f, :set_values)
    s.v_wind_gnd_idx = parameter_index(integrator.f, :v_wind_gnd)
    s.v_wind_idx = parameter_index(integrator.f, :v_wind)
    s.stiffness_factor_idx = parameter_index(integrator.f, :stiffness_factor)
    s.get_pos = getu(integrator.sol, simple_sys.pos[:,:])
    s.get_steering_pos = getu(integrator.sol, simple_sys.steering_pos)
    s.get_line_acc = getu(integrator.sol, simple_sys.acc[:,s.num_E-2])
    s.get_kite_vel = getu(integrator.sol, simple_sys.vel[:,s.num_A])
    s.get_winch_forces = getu(integrator.sol, simple_sys.force[:,1:3])
    s.get_L_C = getu(integrator.sol, simple_sys.L_C)
    s.get_L_D = getu(integrator.sol, simple_sys.L_D)
    s.get_D_C = getu(integrator.sol, simple_sys.D_C)
    s.get_D_D = getu(integrator.sol, simple_sys.D_D)
    s.get_tether_lengths = getu(integrator.sol, simple_sys.tether_length)
    s.get_tether_speeds = getu(integrator.sol, simple_sys.tether_speed)
    update_pos!(s, integrator)
    return integrator
end

# remove this in favor of using init_sim!
function reset_sim!(s::KPS4_3L; stiffness_factor=0.035)
    clear!(s)
    s.stiffness_factor = stiffness_factor  
    dt = 1/s.set.sample_freq
    # 1. KenCarp4
    # TRBDF2, Rodas4P, Rodas5P, Kvaerno5, KenCarp4, radau, QNDF
    integrator = OrdinaryDiffEq.init(s.prob, KenCarp4(autodiff=false); dt=dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
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
    OrdinaryDiffEq.step!(integrator, dt, true)
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

# ==================== end of getter functions ================================================

# not implemented
function spring_forces(s::KPS4_3L)
    forces = zeros(SimFloat, s.num_A)
    for i in 1:s.set.segments*3
        forces[i] =  s.springs[i].c_spring * (norm(s.pos[i+3] - s.pos[i]) - s.segment_lengths[(i-1)%3+1]) * s.stiffness_factor
        if forces[i] > 4000.0
            println("Tether raptures for segment $i !")
        end
    end
    for i in 1:KITE_SPRINGS_3L
        p1 = s.springs[i+s.set.segments*3].p1  # First point nr.
        p2 = s.springs[i+s.set.segments*3].p2  # Second point nr.
        pos1, pos2 = s.pos[p1], s.pos[p2]
        spring = s.springs[i+s.set.segments*3]
        l_0 = spring.length # Unstressed lengthc_spring
        k = spring.c_spring * s.stiffness_factor       # Spring constant 
        segment = pos1 - pos2
        norm1 = norm(segment)
        k1 = 0.25 * k # compression stiffness kite segments
        if (norm1 - l_0) > 0.0
            spring_force = k *  (norm1 - l_0) 
        else 
            spring_force = k1 *  (norm1 - l_0)
        end
        forces[i+s.set.segments*3] = spring_force
        if norm(s.spring_force) > 4000.0
            println("Bridle brakes for spring $i !")
        end
    end
    forces
end


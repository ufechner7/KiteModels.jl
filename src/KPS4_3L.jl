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
    flap_acc::MVector{2, S} = zeros(MVector{2, S})
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
    flap_angle::MVector{2, S} =      zeros(S, 2)
    "air density at the height of the kite"
    rho::S =               0.0
    "multiplier for the damping of all movement"
    damping_coeff::S =  0.0
    "current masses, depending on the total tether length"
    masses::MVector{P, S}         = zeros(P)
    "vector of the springs, defined as struct"
    springs::MVector{Q, SP}       = zeros(SP, Q)
    "unit spring coefficient"
    c_spring::S = zero(S)
    "unit damping coefficient"
    damping::S = zero(S)
    "vector of the forces, acting on the particles"
    forces::SVector{P, T} = zeros(SVector{P, T})
    "synchronous speed or torque of the motor/ generator"
    set_values::KVec3  = zeros(KVec3)
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
    "Lift of point C"
    L_C::T = zeros(S, 3)
    "Lift of point D"
    L_D::T = zeros(S, 3)
    "Drag of point C"
    D_C::T = zeros(S, 3)
    "Drag of point D"
    D_D::T = zeros(S, 3)
    "Solution of the steady state problem"
    steady_sol::Union{SciMLBase.NonlinearSolution, Nothing} = nothing
    "Simplified system of the mtk model"
    simple_sys::Union{ModelingToolkit.ODESystem, Nothing} = nothing

    # set_values_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    v_wind_gnd_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    v_wind_idx::Union{ModelingToolkit.ParameterIndex, Nothing} = nothing
    prob::Union{OrdinaryDiffEqCore.ODEProblem, Nothing} = nothing
    set_set_values::Union{SymbolicIndexingInterface.MultipleSetters, Nothing} = nothing
    get_pos::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_flap_angle::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_flap_acc::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_kite_vel::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_winch_forces::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_tether_lengths::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_tether_vels::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_L_C::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_L_D::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_D_C::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
    get_D_D::Union{SymbolicIndexingInterface.MultipleGetters, SymbolicIndexingInterface.TimeDependentObservedFunction, Nothing} = nothing
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
    s.flap_angle .= 0.0
    s.vel_kite .= 0.0
    [s.winch_forces[i] .= 0.0 for i in 1:3]
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
    for i in 1:s.num_A
        s.forces[i] .= 0.0
        s.veld[i] .= 0.0
    end
    s.rho = s.set.rho_0
    s.c_spring = s.set.e_tether * (s.set.d_tether/2000.0)^2 * pi
    s.damping = (s.set.damping / s.set.c_spring) * s.c_spring
    init_masses!(s)
    init_springs!(s)
end

# include(joinpath(@__DIR__, "CreatePolars.jl"))
function KPS4_3L(kcu::KCU)
    set = kcu.set
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
    cl_struct = linear_interpolation((alphas, d_flap_angles), cl_matrix; extrapolation_bc = NaN)
    cd_struct = linear_interpolation((alphas, d_flap_angles), cd_matrix; extrapolation_bc = NaN)
    c_te_struct = linear_interpolation((alphas, d_flap_angles), c_te_matrix; extrapolation_bc = NaN)
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
    ss.heading = calc_heading(s)
    ss.course = calc_course(s)
    ss.v_app = norm(s.v_apparent)
    ss.l_tether = s.tether_lengths[3]
    ss.v_reelout = s.reel_out_speeds[3]
    ss.depower = rad2deg(s.flap_angle[1] + s.flap_angle[2])
    ss.steering = rad2deg(s.flap_angle[2] - s.flap_angle[1])
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
    depower = rad2deg(s.flap_angle[1] + s.flap_angle[2])
    steering = rad2deg(s.flap_angle[2] - s.flap_angle[1])
    KiteUtils.SysState{P}(s.t_0, t_sim, 0, 0, orient, elevation, azimuth, s.tether_lengths[3], s.reel_out_speeds[3], forces[3], depower, steering, 
                          heading, course, v_app_norm, s.vel_kite, X, Y, Z, 
                          0, 0, 0, 0, 
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end


function calc_heading(s::KPS4_3L)
    # turn s.e_x by -azimuth around global z-axis and then by elevation around global y-axis
    vec = rotate_in_xz(rotate_in_yx(s.e_x, -calc_azimuth(s)), -calc_elevation(s))
    heading = atan(-vec[2], vec[3])
    return heading
end

"""
    init_sim!(s; damping_coeff=1.0, prn=false, torque_control=true)

Initialises the integrator of the model.

Parameters:
- s:     an instance of a 3 line kite model
- damping_coeff: amount of damping in the first steps
- prn: if set to true, print the detailed solver results
- torque_control: wether or not to use torque control

Returns:
Nothing.
"""
function init_sim!(s::KPS4_3L; damping_coeff=50.0, prn=false, 
                   torque_control=true) # TODO: add sysstate init ability
    clear!(s)
    change_control_mode = s.torque_control != torque_control
    s.torque_control = torque_control

    dt = 1/s.set.sample_freq*2
    tspan   = (0.0, dt) 
    solver = FBDF(autodiff=false) # https://docs.sciml.ai/SciMLBenchmarksOutput/stable/#Results
    s.damping_coeff = damping_coeff
    new_inital_conditions = (s.last_init_elevation != s.set.elevation || s.last_init_tether_length != s.set.l_tether)
    s.set_hash = settings_hash(s.set)
    init_new_model = isnothing(s.prob) || change_control_mode || s.last_set_hash != s.set_hash
    init_new_pos = new_inital_conditions && !isnothing(s.get_pos)

    if init_new_model
        if prn; println("initializing with new model and new pos"); end
        pos, vel = init_pos_vel(s)
        sys, inputs = model!(s, pos, vel)
        (s.simple_sys, _) = structural_simplify(sys, (inputs, []); simplify=true)
        s.prob = ODEProblem(s.simple_sys, nothing, tspan; fully_determined=true)
        s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
        next_step!(s; set_values=zeros(3), dt=1.0) # step to get stable state
        s.u0 = deepcopy(s.integrator.u)
        OrdinaryDiffEqCore.reinit!(s.integrator, s.u0)
    elseif init_new_pos
        if prn; println("initializing with last model and new pos"); end
        pos, vel = init_pos_vel(s)
        pos, vel = convert_pos_vel(s, pos, vel)
        defaults = vcat([vcat([s.simple_sys.pos[j, i] => pos[j, i] for i in 1:s.num_flap_C-1 for j in 1:3]), 
                        vcat([s.simple_sys.pos[j, i] => pos[j, i] for i in s.num_flap_D+1:s.num_A for j in 1:3]),
                        s.simple_sys.tether_length => s.tether_lengths]...)
        s.prob = ODEProblem(s.simple_sys, defaults, tspan)
        OrdinaryDiffEqCore.reinit!(s.integrator, s.prob.u0)
        next_step!(s; set_values=zeros(3), dt=1.0) # step to get stable state
        s.u0 = deepcopy(s.integrator.u)
        OrdinaryDiffEqCore.reinit!(s.integrator, s.u0)
    else
        if prn; println("initializing with last model and last pos"); end
        OrdinaryDiffEqCore.reinit!(s.integrator, s.u0)
        # s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
    end

    s.last_init_elevation = s.set.elevation
    s.last_init_tether_length = s.set.l_tether
    s.last_set_hash = s.set_hash
    update_pos!(s)
    return nothing
end


function next_step!(s::KPS4_3L; set_values=zeros(KVec3), v_wind_gnd=s.set.v_wind, upwind_dir=-pi/2, dt=1/s.set.sample_freq)
    s.iter = 0
    set_v_wind_ground!(s, calc_height(s), v_wind_gnd, upwind_dir)
    if isnothing(s.get_pos)
        s.v_wind_gnd_idx = parameter_index(s.integrator.f, :v_wind_gnd)
        s.v_wind_idx = parameter_index(s.integrator.f, :v_wind)
        s.set_set_values = setu(s.integrator.sol, s.simple_sys.set_values)
        s.get_pos = getu(s.integrator.sol, s.simple_sys.pos[:,:])
        s.get_flap_angle = getu(s.integrator.sol, s.simple_sys.flap_angle)
        s.get_flap_acc = getu(s.integrator.sol, s.simple_sys.flap_acc)
        s.get_kite_vel = getu(s.integrator.sol, s.simple_sys.vel[:,s.num_A])
        s.get_winch_forces = getu(s.integrator.sol, s.simple_sys.force[:,1:3])
        s.get_L_C = getu(s.integrator.sol, s.simple_sys.L_C)
        s.get_L_D = getu(s.integrator.sol, s.simple_sys.L_D)
        s.get_D_C = getu(s.integrator.sol, s.simple_sys.D_C)
        s.get_D_D = getu(s.integrator.sol, s.simple_sys.D_D)
        s.get_tether_lengths = getu(s.integrator.sol, s.simple_sys.tether_length)
        s.get_tether_vels = getu(s.integrator.sol, s.simple_sys.tether_vel)
    end
    s.set_values .= set_values
    s.set_set_values(s.integrator, s.set_values)
    s.integrator.p[s.v_wind_gnd_idx] .= s.v_wind_gnd
    s.integrator.p[s.v_wind_idx] .= s.v_wind
    s.t_0 = s.integrator.t
    OrdinaryDiffEqCore.step!(s.integrator, dt, true)
    if !successful_retcode(s.integrator.sol)
        println("Return code for solution: ", s.integrator.sol.retcode)
    end
    @assert successful_retcode(s.integrator.sol)
    update_pos!(s)
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


function calc_acc_speed(motor::AsyncMachine, tether_vel, norm_, set_speed)
    calc_acceleration(motor, tether_vel, norm_; set_speed, set_torque=nothing, use_brake=true)
end
@register_symbolic calc_acc_speed(motor::AsyncMachine, tether_vel, norm_, set_speed)

function calc_acc_torque(motor::TorqueControlledMachine, tether_vel, norm_, set_torque)
    calc_acceleration(motor, tether_vel, norm_; set_speed=nothing, set_torque, use_brake=false)
end
@register_symbolic calc_acc_torque(motor::TorqueControlledMachine, tether_vel, norm_, set_torque)

function sym_interp(interp::Function, aoa, flap_angle)
    return interp(rad2deg(aoa), rad2deg(flap_angle-aoa))
end
@register_symbolic sym_interp(interp::Function, aoa, flap_angle)


"""
    calc_aero_forces!(s::KPS4_3L, eqs2, force_eqs, force, pos, vel, t, e_x, e_y, e_z, rho)

Calculates the aerodynamic forces acting on the kite particles.

Parameters:
- pos:              vector of the particle positions
- vel:              vector of the particle velocities
- rho:              air density [kg/m^3]

Updates the vector s.forces of the first parameter.
"""
function calc_aero_forces_mtk!(s::KPS4_3L, eqs2, force_eqs, force, pos, vel, t, e_x, e_y, e_z, E_C, rho, v_wind, flap_angle)
    n = s.set.aero_surfaces
    @variables begin
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
        v_cx    ~ (vel[:, s.num_C] ⋅ e_x) * e_x
        v_dx    ~ (vel[:, s.num_D] ⋅ e_x) * e_x
        v_dy    ~ (vel[:, s.num_D] ⋅ e_y) * e_y
        v_dz    ~ (vel[:, s.num_D] ⋅ e_z) * e_z
        v_cy    ~ (vel[:, s.num_C] ⋅ e_y) * e_y
        v_cz    ~ (vel[:, s.num_C] ⋅ e_z) * e_z
        y_lc    ~  norm(pos[:, s.num_C] - 0.5 * (pos[:, s.num_C] + pos[:, s.num_D]))
        y_ld    ~ -norm(pos[:, s.num_D] - 0.5 * (pos[:, s.num_C] + pos[:, s.num_D]))
    ]

    # integrating loop variables, iterating over 2n segments
    @variables begin
        F(t)[1:3, 1:2n]
        e_r(t)[1:3, 1:2n]
        e_te(t)[1:3, 1:2n] # clockwise trailing edge of flap vector
        y_l(t)[1:2n]
        v_kite(t)[1:3, 1:2n]
        v_a(t)[1:3, 1:2n]
        e_drift(t)[1:3, 1:2n]
        v_a_xr(t)[1:3, 1:2n]
        aoa(t)[1:n*2]
        cl_seg(t)[1:n*2]
        cd_seg(t)[1:n*2]
        L_seg(t)[1:3, 1:2n]
        D_seg(t)[1:3, 1:2n]
        F_te_seg(t)[1:3, 1:2n]
        seg_flap_angle(t)[1:2n] # flap angle relative to -e_x, e_z
        ram_force(t)[1:2n]
        te_force(t)[1:2n]
        L_C(t)[1:3]
        L_D(t)[1:3]
        D_C(t)[1:3]
        D_D(t)[1:3]
        F_te_C(t)[1:3] # trailing edge clockwise force
        F_te_D(t)[1:3]
    end
    l_c_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(L_C .~ 0))
    l_d_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(L_D .~ 0))
    d_c_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(D_C .~ 0))
    d_d_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(D_D .~ 0))
    f_te_c_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(F_te_C .~ 0))
    f_te_d_eq = SizedArray{Tuple{3}, Symbolics.Equation}(collect(F_te_D .~ 0))
    kite_length = zero(SimFloat)
    α           = zero(SimFloat)
    α_0         = zero(SimFloat)
    α_middle    = zero(SimFloat)
    dα          = zero(SimFloat)
    α_0         = π/2 - s.set.width/2/s.set.radius
    α_middle    = π/2
    dα          = (α_middle - α_0) / n
    for i in 1:n*2
        if i <= n
            α = α_0 + -dα/2 + i * dα
            kite_length = s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (α - α_0) / (π/2 - α_0) # TODO: kite length gets less with flap turning
        else
            α = pi - (α_0 + -dα/2 + (i-n) * dα)
            kite_length = s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (π - α_0 - α) / (π/2 - α_0)
        end
        seg_flap_height = kite_length * s.set.flap_height
        eqs2 = [
            eqs2
            F[:, i]          ~ E_C + e_y * cos(α) * s.set.radius - e_z * sin(α) * s.set.radius
            e_r[:, i]        ~ (E_C - F[:, i]) / norm(E_C - F[:, i])
            y_l[i]           ~ cos(α) * s.set.radius
            α < π/2 ?
                v_kite[:, i] ~ ((v_cx - v_dx) / (y_lc - y_ld) * (y_l[i] - y_ld) + v_dx) + v_cy + v_cz :
                v_kite[:, i] ~ ((v_cx - v_dx) / (y_lc - y_ld) * (y_l[i] - y_ld) + v_dx) + v_dy + v_dz
            v_a[:, i]         ~ v_wind .- v_kite[:, i]
            e_drift[:, i]    ~ (e_r[:, i] × e_x)
            v_a_xr[:, i]     ~ v_a[:, i] .- (v_a[:, i] ⋅ e_drift[:, i]) .* e_drift[:, i]

            α < s.α_l ?
                seg_flap_angle[i]    ~ flap_angle[1]  :
            α > s.α_r ?
                seg_flap_angle[i]    ~ flap_angle[2] :
                seg_flap_angle[i]    ~ ((flap_angle[2] - flap_angle[1]) / (s.α_r - s.α_l) * (α - s.α_l) + (flap_angle[1]))

            aoa[i]      ~ -asin((v_a_xr[:, i] / norm(v_a_xr[:, i])) ⋅ e_r[:, i]) + deg2rad(s.set.alpha_zero)
            cl_seg[i]   ~ sym_interp(s.cl_interp, aoa[i], seg_flap_angle[i])
            cd_seg[i]   ~ sym_interp(s.cd_interp, aoa[i], seg_flap_angle[i])

            L_seg[:, i] ~ 0.5 * rho * (norm(v_a_xr[:, i]))^2 * s.set.radius * dα * kite_length * cl_seg[i] * 
                                ((v_a_xr[:, i] × e_drift[:, i]) / norm(v_a_xr[:, i] × e_drift[:, i]))
            D_seg[:, i] ~ 0.5 * rho * norm(v_a_xr[:, i]) * s.set.radius * dα * kite_length * cd_seg[i] *
                                v_a_xr[:, i]


            e_te[:, i] ~ e_x * sin(seg_flap_angle[i]) + e_r[:, i] * cos(seg_flap_angle[i])
            ram_force[i] ~ smooth_sign(deg2rad(s.set.alpha_zero) - seg_flap_angle[i]) *
                        rho * norm(v_a[:, i])^2 * seg_flap_height * s.set.radius * dα * (seg_flap_height/2) / (kite_length/4)
            te_force[i] ~ 0.5 * rho * (norm(v_a_xr[:, i]))^2 * s.set.radius * dα * kite_length * 
                                sym_interp(s.c_te_interp, aoa[i], seg_flap_angle[i])
            F_te_seg[:, i] ~ (ram_force[i] + te_force[i]) * e_te[:, i]
        ]

        # TODO: correct for extra torque in wingtips (add to c substract from d)
        # TODO: use SymbolicNumericIntegration.jl
        if i <= n
            [l_c_eq[j] = (L_C[j] ~ l_c_eq[j].rhs + L_seg[j, i]) for j in 1:3]
            [d_c_eq[j] = (D_C[j] ~ d_c_eq[j].rhs + D_seg[j, i]) for j in 1:3]
            [f_te_c_eq[j] = (F_te_C[j] ~ f_te_c_eq[j].rhs + F_te_seg[j, i]) for j in 1:3]
        else 
            [l_d_eq[j] = (L_D[j] ~ l_d_eq[j].rhs + L_seg[j, i]) for j in 1:3]
            [d_d_eq[j] = (D_D[j] ~ d_d_eq[j].rhs + D_seg[j, i]) for j in 1:3]
            [f_te_d_eq[j] = (F_te_D[j] ~ f_te_d_eq[j].rhs + F_te_seg[j, i]) for j in 1:3]
        end
    end

    
    eqs2 = [
        eqs2
        l_c_eq
        d_c_eq
        l_d_eq
        d_d_eq
        f_te_c_eq
        f_te_d_eq
    ]

    # longtitudinal force
    # F_inside_flap = P * A
    # F_inside_flap = rho * norm(v)^2 * flap_height * width
    # F_inside_flap = rho * norm(v)^2 * flap_height * radius * dα
    # F_trailing_edge = -F_inside_flap * (flap_height/2) / (kite_length/4) if flap_angle > 0 clockwise force
    # F_trailing_edge = F_inside_flap * (flap_height/2) / (kite_length/4) if flap_angle < 0 clockwise force
    # dF_te_dα = rho * norm(v)^2 * flap_height * radius
    # flap_height = height_middle * kite_length / middle_length
    
    # TODO: check if the right forces are added
    force_eqs[:,s.num_C]   .= (force[:,s.num_C]   .~ L_C + D_C) 
    force_eqs[:,s.num_D]   .= (force[:,s.num_D]   .~ L_D + D_D)
    force_eqs[:,s.num_flap_C] .= (force[:,s.num_flap_C] .~ F_te_C)
    force_eqs[:,s.num_flap_D] .= (force[:,s.num_flap_D] .~ F_te_D)
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
    damping, rho, i, l_0, k, c, segment, rel_vel, av_vel, norm1, unit_vector, k1, k2, c1, c2, spring_vel, perp_vel,
            spring_force, v_apparent, v_wind_tether, area, v_app_perp, half_drag_force)
    d_tether = s.set.d_tether/1000.0
    eqs2 = [
        eqs2
        i <= s.set.segments*3 ? l_0 ~ length[(i-1) % 3 + 1] : l_0 ~ s.springs[i].length # Unstressed length
        i <= s.set.segments*3 ? k   ~ c_spring[(i-1) % 3 + 1] :
                                k   ~ s.springs[i].c_spring        # Spring constant
        i <= s.set.segments*3 ? c   ~ damping[(i-1) % 3 + 1] : c ~ s.springs[i].damping # Damping coefficient    
        segment     .~ pos1 - pos2 # TODO: all segments have same length and tension
        rel_vel     .~ vel1 - vel2
        av_vel      .~ 0.5 * (vel1 + vel2)
        norm1        ~ norm(segment)
        unit_vector .~ segment / norm1
        k1           ~ 1.0 * k # compression stiffness kite segments
        k2           ~ 0.1 * k  # compression stiffness tether segments
        c1           ~ 6.0 * c  # damping kite segments
        c2           ~ 0.05 * c  # damping perpendicular
        spring_vel   ~ rel_vel ⋅ unit_vector
        perp_vel    .~ rel_vel .- spring_vel * unit_vector
    ]

    if i >= Base.length(s.springs) - KITE_SPRINGS_3L + 1  # kite springs
        for j in 1:3
            eqs2 = [
                eqs2
                spring_force[j] ~ 
                    (k  * (l_0 - norm1) - c1 * spring_vel) * unit_vector[j] * (1 + smooth_sign(norm1 - l_0)) / 2 +
                    (k1 * (l_0 - norm1) -  c * spring_vel) * unit_vector[j] * (1 - smooth_sign(norm1 - l_0)) / 2
            ]
        end
    else
        for j in 1:3
            eqs2 = [
                eqs2
                spring_force[j] ~
                    ((k  * (l_0 - norm1) - c * spring_vel) * unit_vector[j] - c2 * perp_vel[j]) * (1 + smooth_sign(norm1 - l_0)) / 2 +
                    ((k2 * (l_0 - norm1) - c * spring_vel) * unit_vector[j] - c2 * perp_vel[j]) * (1 - smooth_sign(norm1 - l_0)) / 2
            ]
        end
    end
    eqs2 = [
        eqs2
        v_apparent       ~ v_wind_tether - av_vel
        i >= s.num_flap_C ?
            area             ~ norm1 * d_tether * 10 : # 10 is the number of parallel lines in the bridle system
            area             ~ norm1 * d_tether * (1 + (i%3 == 0)) # double area for middle tether
        v_app_perp       ~ v_apparent - (v_apparent ⋅ unit_vector) * unit_vector
        half_drag_force .~ (0.25 * rho * s.set.cd_tether * norm(v_app_perp) * area) .* v_app_perp
    ]

    for j in 1:3
        force_eqs[j, s.springs[i].p1] = 
            (force[j, s.springs[i].p1] ~ force_eqs[j, s.springs[i].p1].rhs + (half_drag_force[j] + spring_force[j]))
        force_eqs[j, s.springs[i].p2] = 
            (force[j, s.springs[i].p2] ~ force_eqs[j, s.springs[i].p2].rhs + (half_drag_force[j] - spring_force[j]))
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
@inline function inner_loop_mtk!(s::KPS4_3L, eqs2, force_eqs, t, force, pos, vel, length, c_spring, damping, v_wind_gnd)
    @variables begin
        height(t)[eachindex(s.springs)]
        rho(t)[eachindex(s.springs)]
        v_wind_tether(t)[1:3, eachindex(s.springs)]

        l_0(t)[eachindex(s.springs)]
        k(t)[eachindex(s.springs)]
        c(t)[eachindex(s.springs)]
        segment(t)[1:3, eachindex(s.springs)]
        rel_vel(t)[1:3, eachindex(s.springs)]
        av_vel(t)[1:3, eachindex(s.springs)] 
        norm1(t)[eachindex(s.springs)]
        unit_vector(t)[1:3, eachindex(s.springs)]
        k1(t)[eachindex(s.springs)]
        k2(t)[eachindex(s.springs)]
        c1(t)[eachindex(s.springs)]
        c2(t)[eachindex(s.springs)]
        spring_vel(t)[eachindex(s.springs)]
        perp_vel(t)[1:3, eachindex(s.springs)]
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
            height[i]           ~ max(0.0, 0.5 * (pos[:, p1][3] + pos[:, p2][3]))
            rho[i]              ~ calc_rho(s.am, height[i])
            v_wind_tether[:, i] ~ calc_wind_factor(s.am, height[i]) * v_wind_gnd
        ]

        eqs2, force_eqs = calc_particle_forces_mtk!(s, eqs2, force_eqs, force, pos[:, p1], pos[:, p2], vel[:, p1], 
                          vel[:, p2], length, c_spring, damping, rho[i], i, l_0[i], k[i], c[i], segment[:, i], 
                          rel_vel[:, i], av_vel[:, i], norm1[i], unit_vector[:, i], k1[i], k2[i], c1[i], c2[i], spring_vel[i],
                          perp_vel[:, i], spring_force[:, i], v_apparent[:, i], v_wind_tether[:, i], area[i], v_app_perp[:, i], 
                          half_drag_force[:, i])
    end

    return eqs2, force_eqs
end

function update_pos!(s)
    pos = s.get_pos(s.integrator)
    s.flap_angle       .= s.get_flap_angle(s.integrator)
    s.flap_acc         .= s.get_flap_acc(s.integrator)
    [s.pos[i]          .= pos[:, i] for i in 1:s.num_A]
    s.vel_kite         .= s.get_kite_vel(s.integrator)
    winch_forces        = s.get_winch_forces(s.integrator)
    [s.winch_forces[i] .= (winch_forces[:, i]) for i in 1:3]
    s.tether_lengths   .= s.get_tether_lengths(s.integrator)
    s.reel_out_speeds  .= s.get_tether_vels(s.integrator)
    s.L_C               = s.get_L_C(s.integrator)
    s.L_D               = s.get_L_D(s.integrator)
    s.D_C               = s.get_D_C(s.integrator)
    s.D_D               = s.get_D_D(s.integrator)
    calc_kite_ref_frame!(s, s.pos[s.num_E], s.pos[s.num_C], s.pos[s.num_D])
    # @assert all(abs.(s.flap_angle) .<= deg2rad(90))
    nothing
end

function convert_pos_vel(s::KPS4_3L, pos_, vel_)
    pos = Array{Union{Nothing, Float64}}(nothing, 3, s.num_A)
    vel = Array{Union{Nothing, Float64}}(nothing, 3, s.num_A)
    [pos[:,i] .= pos_[i] for i in 1:s.num_flap_C-1]
    [vel[:,i] .= vel_[i] for i in 1:s.num_flap_C-1]
    [pos[:,i] .= pos_[i] for i in s.num_flap_D+1:s.num_A]
    [vel[:,i] .= vel_[i] for i in s.num_flap_D+1:s.num_A]
    return pos, vel
end

function model!(s::KPS4_3L, pos_, vel_)
    pos_, vel_ = convert_pos_vel(s, pos_, vel_)
    if s.torque_control
        [s.motors[i] = TorqueControlledMachine(s.set) for i in 1:3]
    else
        [s.motors[i] = AsyncMachine(s.set) for i in 1:3]
    end
    @parameters begin
        v_wind_gnd[1:3] = s.v_wind_gnd
        v_wind[1:3] = s.v_wind
    end
    @variables begin
        set_values(t)[1:3] = s.set_values
        pos(t)[1:3, 1:s.num_A] = pos_
        vel(t)[1:3, 1:s.num_A] = vel_
        acc(t)[1:3, 1:s.num_A]
        flap_angle(t)[1:2]   = zeros(2) # angle
        flap_vel(t)[1:2]     = zeros(2) # angular vel
        flap_acc(t)[1:2]                # angular acc
        tether_length(t)[1:3]  = s.tether_lengths
        tether_vel(t)[1:3] = zeros(3)
        segment_length(t)[1:3]
        mass_tether_particle(t)[1:3]
        damping(t)[1:3]
        damping_coeff(t)
        c_spring(t)[1:3]
        P_c(t)[1:3]
        E_C(t)[1:3]
        e_x(t)[1:3]
        e_y(t)[1:3]
        e_z(t)[1:3]
        e_r_C(t)[1:3]
        e_r_D(t)[1:3]
        e_te_C(t)[1:3]
        e_te_D(t)[1:3]
        force(t)[1:3, 1:s.num_A]
        rho_kite(t)
    end
    # Collect the arrays into variables
    pos = collect(pos)
    vel = collect(vel)
    acc = collect(acc)

    eqs1 = []
    mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2

    [eqs1 = vcat(eqs1, D.(pos[:, i]) .~ 0.0) for i in 1:3]
    [eqs1 = vcat(eqs1, D.(pos[:, i]) .~ vel[:, i]) for i in 4:s.num_flap_C-1]
    eqs1 = [eqs1; D(flap_angle)   ~ flap_vel]
    [eqs1 = vcat(eqs1, D.(pos[:, i]) .~ vel[:, i]) for i in s.num_E:s.num_A]
    [eqs1 = vcat(eqs1, D.(vel[:, i]) .~ 0.0) for i in 1:3]
    [eqs1 = vcat(eqs1, D.(vel[:, i]) .~ acc[:, i]) for i in 4:s.num_flap_C-1]
    eqs1 = [eqs1; D(flap_vel)   ~ flap_acc]
    [eqs1 = vcat(eqs1, D.(vel[:, i]) .~ acc[:, i]) for i in s.num_E:s.num_A]

    eqs1 = vcat(eqs1, D.(tether_length) .~ tether_vel)
    if s.torque_control
        eqs1 = vcat(eqs1, D.(tether_vel) .~ [calc_acc_torque(s.motors[i], tether_vel[i], norm(force[:, (i-1) % 3 + 1]),
                                                               set_values[i]) for i in 1:3])
    else
        eqs1 = vcat(eqs1, D.(tether_vel) .~ [calc_acc_speed(s.motors[i], tether_vel[i], norm(force[:,(i-1) % 3 + 1]), 
                                                               set_values[i]) for i in 1:3])
    end

    # Compute the masses and forces
    force_eqs = SizedArray{Tuple{3, s.num_A}, Symbolics.Equation}(undef)
    force_eqs[:, :] .= (force[:, :] .~ 0)

    flap_length = s.kite_length_C/4
    eqs2 = [
        pos[:, s.num_flap_C]    ~ pos[:, s.num_C] - e_x * flap_length * cos(flap_angle[1]) + e_r_C * flap_length * sin(flap_angle[1])
        pos[:, s.num_flap_D]    ~ pos[:, s.num_D] - e_x * flap_length * cos(flap_angle[2]) + e_r_D * flap_length * sin(flap_angle[2])
        vel[:, s.num_flap_C]    ~ vel[:, s.num_C] - e_x * flap_length * cos(flap_vel[1]) + e_r_C * flap_length * sin(flap_vel[1])
        vel[:, s.num_flap_D]    ~ vel[:, s.num_D] - e_x * flap_length * cos(flap_vel[2]) + e_r_D * flap_length * sin(flap_vel[2])
        acc[:, s.num_flap_C]    ~ acc[:, s.num_C] - e_x * flap_length * cos(flap_acc[1]) + e_r_C * flap_length * sin(flap_acc[1])
        acc[:, s.num_flap_D]    ~ acc[:, s.num_D] - e_x * flap_length * cos(flap_acc[2]) + e_r_D * flap_length * sin(flap_acc[2])
        segment_length          ~ tether_length  ./ s.set.segments
        mass_tether_particle    ~ mass_per_meter .* segment_length
        damping                 ~ [s.damping / segment_length[1], s.damping / segment_length[2], s.damping*2 / segment_length[3]]
        c_spring                ~ [s.c_spring / segment_length[1], s.c_spring / segment_length[2], s.c_spring*2 / segment_length[3]]
        P_c     ~ 0.5 * (pos[:, s.num_C] + pos[:, s.num_D])
        e_y     ~ (pos[:, s.num_C] - pos[:, s.num_D]) / norm(pos[:, s.num_C] - pos[:, s.num_D])
        e_z     ~ (pos[:, s.num_E] - P_c) / norm(pos[:, s.num_E] - P_c)
        e_x     ~ cross(e_y, e_z)
        e_r_C   ~ (E_C - pos[:, s.num_C]) / norm(E_C - pos[:, s.num_C])
        e_r_D   ~ (E_C - pos[:, s.num_D]) / norm(E_C - pos[:, s.num_D])
        e_te_C ~ e_x * sin(flap_angle[1]) + e_r_C * cos(flap_angle[1])
        e_te_D ~ e_x * sin(flap_angle[2]) + e_r_D * cos(flap_angle[2])
        # E_C is the center of the circle shape of the front view of the kite
        E_C     ~ pos[:, s.num_E] + e_z * (-s.set.bridle_center_distance + s.set.radius) 
        rho_kite ~ calc_rho(s.am, pos[3,s.num_A])
        damping_coeff ~ max(1.0 - t/2, 0.0) * s.damping_coeff
    ]

    eqs2, force_eqs = calc_aero_forces_mtk!(s, eqs2, force_eqs, force, pos, vel, t, e_x, e_y, e_z, E_C, rho_kite, v_wind, flap_angle)
    eqs2, force_eqs = inner_loop_mtk!(s, eqs2, force_eqs, t, force, pos, vel, segment_length, c_spring, damping, 
                                      v_wind_gnd)
    
    for i in 1:3
        eqs2 = vcat(eqs2, vcat(force_eqs[:, i]))
        eqs2 = vcat(eqs2, acc[:, i] .~ 0)
    end
    for i in 4:s.num_flap_C-1
        eqs2 = vcat(eqs2, vcat(force_eqs[:, i]))
        eqs2 = vcat(eqs2, acc[:, i] .~ [0.0; 0.0; -G_EARTH] .+ (force[:, i] ./ mass_tether_particle[(i-1)%3+1]) .- damping_coeff * vel[:, i])
    end

    # torque = I * flap_acc
    # flap_acc = torque / (1/3 * (kite_mass/8) * kite_length_c^2)
    # torque = force[:, i] * kite_length_c
    # flap_acc = force[:, i] * kite_length_c / (1/3 * (kite_mass/8) * kite_length_c^2)

    # 1. add all flap + spring + drag forces to flap_C point
    # 2. remove forces not in e_flap_c direction
    # 3. substract forces on point flap_C from point C
    # 4. calculate acceleration from force flap c in e_flap_c direction
    [force_eqs[j, s.num_C] = force[j, s.num_C] ~ force_eqs[j, s.num_C].rhs - force_eqs[j, s.num_flap_C].rhs for j in 1:3]
    [force_eqs[j, s.num_D] = force[j, s.num_D] ~ force_eqs[j, s.num_D].rhs - force_eqs[j, s.num_flap_D].rhs for j in 1:3]
    eqs2 = [
        eqs2
        vcat(force_eqs[:, s.num_flap_C])
        vcat(force_eqs[:, s.num_flap_D])
        flap_acc[1] ~ ((force[:, s.num_flap_C] + [0.0, 0.0, -G_EARTH]) ⋅ e_te_C - s.damping * 0.75 * flap_vel[1]) * # TODO: add turning drag instead of damping
                    flap_length / (1/3 * (s.set.mass/8) * flap_length^2) - (damping_coeff*200) * flap_vel[1]
        flap_acc[2] ~ ((force[:, s.num_flap_D] + [0.0, 0.0, -G_EARTH]) ⋅ e_te_D - s.damping * 0.75 * flap_vel[2]) * 
                    flap_length / (1/3 * (s.set.mass/8) * flap_length^2) - (damping_coeff*200) * flap_vel[2]
    ]

    for i in s.num_E:s.num_A
        eqs2 = vcat(eqs2, vcat(force_eqs[:, i]))
        eqs2 = vcat(eqs2, acc[:, i] .~ [0.0; 0.0; -G_EARTH] .+ (force[:, i] ./ s.masses[i]) .- (damping_coeff) * vel[:, i])
    end

    eqs = vcat(eqs1, eqs2)

    @named sys = ODESystem(Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs))), t)
    return sys, collect(set_values)
end


# ====================== helper functions ====================================

function smooth_sign(x; ϵ = 0.01)
    return x / √(x^2 + ϵ^2)
end
@register_symbolic smooth_sign(x)

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
    distance = 10
    for i in 1:rows
        for j in 1:cols
            if isnan(matrix[i, j])
                neighbors = []
                for d in 1:distance
                    found = false
                    if i-d >= 1 && !isnan(matrix[i-d, j]);
                        push!(neighbors, matrix[i-1, j])
                        found = true
                    end
                    if i+d <= rows && !isnan(matrix[i+d, j])
                        push!(neighbors, matrix[i+1, j])
                        found = true
                    end
                    if j-d >= 1 && !isnan(matrix[i, j-d])
                        push!(neighbors, matrix[i, j-1])
                        found = true
                    end
                    if j+d <= cols && !isnan(matrix[i, j+d])
                        push!(neighbors, matrix[i, j+1])
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
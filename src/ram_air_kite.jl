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

const MeasureFloat = Float32

@with_kw mutable struct Measurement
    set_values::MVector{3, MeasureFloat}    = [-50., -1., -1.]
    tether_length::MVector{3, MeasureFloat} = [51., 51., 51.]
    tether_vel::MVector{3, MeasureFloat}    = zeros(MeasureFloat, 3)
    tether_acc::MVector{3, MeasureFloat}    = zeros(MeasureFloat, 3)
    tether_force::MVector{3, MeasureFloat}  = [540., 3., 3.]
    "elevation and azimuth in spherical coordinate system with columns (left, right) and rows (elevation, azimuth)"
    sphere_pos::Matrix{MeasureFloat}            = deg2rad.([89.0 89.0; 1.0 -1.0])
    sphere_vel::Matrix{MeasureFloat}            = zeros(MeasureFloat, 2, 2)
    sphere_acc::Matrix{MeasureFloat}            = zeros(MeasureFloat, 2, 2)
    "positive azimuth wind direction in right-handed ENU frame relative to east / x-axis"
    wind_dir_gnd::MeasureFloat                  = zero(MeasureFloat)
end

function Base.getproperty(m::Measurement, val::Symbol)
    if val === :elevation
        sphere_pos = getfield(m, :sphere_pos)
        return 0.5(sphere_pos[1, 1] + sphere_pos[1, 2])
    elseif val === :azimuth
        sphere_pos = getfield(m, :sphere_pos)
        return 0.5(sphere_pos[2, 1] + sphere_pos[2, 2])
    else
        return getfield(m, val)
    end
end

@enum SegmentType begin
    POWER
    STEERING
    BRIDLE
end

@enum DynamicsType begin
    DYNAMIC
    STATIC
    KITE
    WINCH
end

"""
A normal freely moving tether point
"""
mutable struct Point
    idx::Int16
    pos_b::KVec3 # pos relative to kite COM in body frame
    pos_w::KVec3 # pos in world frame
    type::DynamicsType
end
function Point(idx, pos_b, type)
    Point(idx, pos_b, copy(pos_b), type)
end

"""
Set of bridle lines that share the same twist angle and trailing edge angle
"""
struct KitePointGroup
    idx::Int16
    points::Vector{Int16}
    fixed_index::Int16 # point which the group rotates around under kite deformation
    chord::KVec3 # chord vector in body frame which the group rotates around under kite deformation
    y_airf::KVec3 # spanwise vector in local panel frame which the group rotates around under kite deformation
    type::DynamicsType
end

"""
A segment from one point index to another point index
"""
mutable struct Segment
    idx::Int16
    points::Tuple{Int16, Int16}
    type::SegmentType
    l0::SimFloat
    diameter::SimFloat
end
function Segment(idx, points, type)
    Segment(idx, points, type, zero(SimFloat), zero(SimFloat))
end
function Segment(idx, points, type, l0)
    Segment(idx, points, type, l0, zero(SimFloat))
end

"""
A pulley described by two segments with the common point of the segments being the pulley
"""
mutable struct Pulley
    idx::Int16
    segments::Tuple{Int16, Int16}
    type::DynamicsType
    sum_length::SimFloat
    function Pulley(idx, segments, type)
        new(idx, segments, type, zero(SimFloat))
    end
end

"""
A set of segments making a flexible tether. The winch point should only be part of one segment.
"""
struct Tether
    idx::Int16
    segments::Vector{Int16}
    winch_point::Int16
end

"""
A set of tethers or just one tether connected to a winch
"""
mutable struct Winch
    idx::Int16
    model::AbstractWinchModel
    tethers::Vector{Int16}
    tether_length::Float64
    function Winch(idx, model, tethers)
        new(idx, model, tethers, zero(Float64))
    end
end

struct PointMassSystem
    points::Vector{Point}
    groups::Vector{KitePointGroup}
    segments::Vector{Segment}
    pulleys::Vector{Pulley}
    tethers::Vector{Tether}
    winches::Vector{Winch}
    function PointMassSystem(points, groups, segments, pulleys, tethers, winches)
        for (i, point) in enumerate(points)
            @assert point.idx == i
        end
        for (i, group) in enumerate(groups)
            @assert group.idx == i
        end
        for (i, segment) in enumerate(segments)
            @assert segment.idx == i
        end
        for (i, pulley) in enumerate(pulleys)
            @assert pulley.idx == i
        end
        for (i, tether) in enumerate(tethers)
            @assert tether.idx == i
        end
        for (i, winch) in enumerate(winches)
            @assert winch.idx == i
        end
        new(points, groups, segments, pulleys, tethers, winches)
    end
end

"""
    mutable struct RamAirKite{S, T, P, Q, SP} <: AbstractKiteModel

State of the kite power system, using a quaternion kite model and three steering lines to the ground. Parameters:
- S: Scalar type, e.g. SimFloat
  In the documentation mentioned as Any, but when used in this module it is always SimFloat and not Any.
- V: Vector type, e.g. KVec3
- P: number of tether points of the system, 3segments+3
Normally a user of this package will not have to access any of the members of this type directly,
use the input and output functions instead.

$(TYPEDFIELDS)
"""
@with_kw mutable struct RamAirKite{S, V, P} <: AbstractKiteModel
    "Reference to the settings struct"
    set::Settings
    "Reference to the geometric wing model"
    wing::VortexStepMethod.RamAirWing
    "Reference to the aerodynamic wing model"
    aero::VortexStepMethod.BodyAerodynamics
    "Reference to the VSM aerodynamics solver"
    vsm_solver::VortexStepMethod.Solver
    "Reference to the point mass system with points, segments, pulleys and tethers"
    point_system::PointMassSystem
    "Reference to the atmospheric model as implemented in the package AtmosphericModels"
    am::AtmosphericModel = AtmosphericModel()
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
    "Simplified system of the mtk model"
    sys::Union{ModelingToolkit.ODESystem, Nothing} = nothing
    "Velocity of the kite"
    vel_kite::V =           zeros(S, 3)
    "Inertia around kite x y and z axis of the body frame"
    I_b::V = zeros(S, 3)
    "Initialization values for kite state"
    u0map::Union{Vector{Pair{Num, S}}, Nothing} = nothing
    "Initialization values for kite parameters"
    p0map::Union{Vector{Pair{Num, S}}, Nothing} = nothing
    "X coordinate on normalized 2d foil of bridle attachments"
    bridle_fracs::V = [0.088, 0.31, 0.58, 0.93]
    crease_frac::S = 0.82
    "The top bridle points that are not on the kite, in CAD frame"
    top_bridle_points::Vector{V} = [[0.290199, 0.784697, -2.61305], [0.392683, 0.785271, -2.61201], [0.498202, 0.786175, -2.62148], [0.535543, 0.786175, -2.62148]]
    "Tether diameter of tethers in bridle system [mm]"
    bridle_tether_diameter::SimFloat = 2.
    "Tether diameter of the power tethers [mm]"
    power_tether_diameter::SimFloat = 2.
    "Tether diameter of the steering tethers [mm]"
    steering_tether_diameter::SimFloat = 1.
    "Number of solve! calls"
    iter::Int64 = 0

    unknowns_vec::Vector{SimFloat} = zeros(SimFloat, 3)

    set_set_values::Function       = () -> nothing
    set_measure::Function          = () -> nothing
    set_vsm::Function              = () -> nothing
    set_unknowns::Function             = () -> nothing

    get_state::Function            = () -> nothing
    get_y::Function                = () -> nothing

    prob::Union{OrdinaryDiffEqCore.ODEProblem, Nothing} = nothing
    init_prob::Union{OrdinaryDiffEqCore.ODEProblem, Nothing} = nothing
    integrator::Union{OrdinaryDiffEqCore.ODEIntegrator, Sundials.CVODEIntegrator, Nothing} = nothing
    init_integrator::Union{OrdinaryDiffEqCore.ODEIntegrator, Sundials.CVODEIntegrator, Nothing} = nothing
end

function RamAirKite(set::Settings, aero::BodyAerodynamics, vsm_solver::VortexStepMethod.Solver, point_system::PointMassSystem)
    length(aero.wings) > 1 && throw(ArgumentError("Just one wing allowed in BodyAerodynamics object"))
    wing = aero.wings[1]
    if set.winch_model == "TorqueControlledMachine"
        s = RamAirKite{SimFloat, Vector{SimFloat}, 3*(set.segments + 1)}(
            ; set, wing, aero, vsm_solver, point_system
            )
        s.torque_control = true
    else
        s = RamAirKite{SimFloat, Vector{SimFloat}, 3*(set.segments + 1)}(
            ; set, wing, aero, vsm_solver, point_system
            )
        s.torque_control = false
    end
    return s
end

function update_sys_state!(ss::SysState, s::RamAirKite, zoom=1.0)
    ss.time = s.t_0
    pos, acc, Q_b_w, elevation, azimuth, course, heading, e_x, tether_vel, twist, kite_vel = s.get_state(s.integrator)
    P = length(s.point_system.points)
    for i in 1:P
        ss.X[i] = pos[1, i] * zoom
        ss.Y[i] = pos[2, i] * zoom
        ss.Z[i] = pos[3, i] * zoom
    end
    # TODO
    # ss.kite_acc      .= kite_acc(s)
    # ss.left_tether_vel = tether_vel[1]
    # ss.right_tether_vel = tether_vel[2]
    ss.acc = norm(acc)
    ss.orient .= Q_b_w
    ss.elevation = elevation
    ss.azimuth = azimuth
    ss.force = zero(SimFloat)
    ss.heading = heading
    ss.course = course
    ss.l_tether = s.tether_lengths[3]
    ss.v_reelout = mean(tether_vel)
    ss.depower = rad2deg(mean(twist))
    ss.steering = rad2deg(twist[end] - twist[1])
    ss.vel_kite .= kite_vel
    nothing
end

function SysState(s::RamAirKite, zoom=1.0) # TODO: add left and right lines, stop using getters and setters
    isnothing(s.integrator) && throw(ArgumentError("run init_sim!(s) first"))
    pos, acc, Q_b_w, elevation, azimuth, course, heading, e_x, tether_vel, twist, kite_vel = s.get_state(s.integrator)
    P = length(s.point_system.points)
    X = zeros(MVector{P, MyFloat})
    Y = zeros(MVector{P, MyFloat})
    Z = zeros(MVector{P, MyFloat})
    for i in 1:P
        X[i] = pos[1, i] * zoom
        Y[i] = pos[2, i] * zoom
        Z[i] = pos[3, i] * zoom
    end
    
    orient = MVector{4, Float32}(Q_b_w) # TODO: add Q_b_w
    # forces = s.get_tether_force() # TODO: add tether force
    forces = zeros(3)
    t_sim = 0
    depower = rad2deg(mean(twist))
    steering = rad2deg(twist[end] - twist[1])
    ss = SysState{P}()
    ss.time = s.t_0
    ss.t_sim = t_sim
    ss.orient .= orient
    ss.elevation = elevation
    ss.azimuth = azimuth
    ss.l_tether = s.tether_lengths[3]
    ss.v_reelout = tether_vel[3]
    ss.force = forces[3]
    ss.depower = depower
    ss.steering = steering
    ss.heading = heading
    ss.course = course
    ss.vel_kite .= kite_vel
    ss.X = X
    ss.Y = Y
    ss.Z = Z
    ss
end

"""
    init_sim!(s::RamAirKite; prn=true) -> Nothing

Initialize a complete kite power system model from scratch.

This function performs the following operations:
1. Converts measurement data to quaternion orientation
2. Initializes the point mass system representing the kite and tethers
3. Creates the symbolic MTK system with all equations
4. Simplifies the system equations
5. Creates an ODEProblem
6. Serializes the problem to disk for future reuse
7. Calls `reinit!` to set up the integrator

This is a computationally expensive operation and should only be called when the model
structure changes. For normal simulations, prefer calling `reinit!` directly if a serialized
problem already exists.

# Arguments
- `s::RamAirKite`: The kite power system state object
- `prn::Bool=true`: Whether to print progress information

# Returns
- `Nothing`
"""
function init_sim!(s::RamAirKite, measure::Measurement; prn=true, precompile=false, lin_sys=false)
    function init(s, measure)
        init_Q_b_w, R_b_w = measure_to_q(measure)
        init_kite_pos = init!(s.point_system, s.set, R_b_w)

        init_va = R_b_w' * [s.set.v_wind, 0., 0.]
        
        sys, defaults, guesses, inputs = create_sys!(s, s.point_system, measure; init_Q_b_w, init_kite_pos, init_va, lin_sys)
        prn && @info "Simplifying the system"
        @time sys, _ = structural_simplify(sys, (inputs, []); additional_passes=[ModelingToolkit.IfLifting])
        s.sys = sys
        dt = SimFloat(1/s.set.sample_freq)
        if prn
            @info "Creating ODEProblem"
            @time s.prob = ODEProblem(s.sys, defaults, (0.0, dt); guesses)
        else
            s.prob = ODEProblem(s.sys, defaults, (0.0, dt); guesses)
        end
        serialize(prob_path, s.prob)
        s.integrator = nothing
    end
    prob_path = joinpath(KiteUtils.get_data_path(), get_prob_name(s.set; precompile))
    if !ispath(prob_path)
        init(s, measure)
    end
    try
        reinit!(s, measure)
    catch e
        rm(prob_path)
        @info "Rebuilding the system. This can take some minutes..."
        init(s, measure)
        reinit!(s, measure)
    end
    return nothing
end

function init_sim!(::RamAirKite; prn=true)
    throw(ArgumentError("Use the function init_sim!(s::RamAirKite, measure::Measurement) instead."))
end

"""
    reinit!(s::RamAirKite; prn=true) -> Nothing

Reinitialize an existing kite power system model with new state values.

This function performs the following operations:
1. If no integrator exists yet:
   - Loads a serialized ODEProblem from disk
   - Initializes a new ODE integrator 
   - Generates getter/setter functions for the system
2. Converts measurement data to quaternion orientation
3. Initializes the point mass system with new positions
4. Sets initial values for all state variables
5. Reinitializes the ODE integrator
6. Updates the linearized aerodynamic model

This is more efficient than `init!` as it reuses the existing model structure
and only updates the state variables to match the current `measure`.

# Arguments
- `s::RamAirKite`: The kite power system state object
- `prn::Bool=true`: Whether to print progress information

# Returns
- `Nothing`

# Throws
- `ArgumentError`: If no serialized problem exists (run `init_sim!` first)
"""
function reinit!(s::RamAirKite, measure::Measurement; prn=true, reload=true)
    isnothing(s.point_system) && throw(ArgumentError("PointMassSystem not defined"))

    init_Q_b_w, R_b_w = measure_to_q(measure)
    init_kite_pos = init!(s.point_system, s.set, R_b_w)
    
    if isnothing(s.prob)
        prob_path = joinpath(KiteUtils.get_data_path(), get_prob_name(s.set))
        !ispath(prob_path) && throw(ArgumentError("$prob_path not found. Run init_sim!(s::RamAirKite) first."))
        try
            s.prob = deserialize(prob_path)
        catch e
            @warn "Failure to deserialize $prob_path !"
            throw(e)
        end
    end
    if isnothing(s.integrator) || !successful_retcode(s.integrator.sol) || reload
        t = @elapsed begin
            dt = SimFloat(1/s.set.sample_freq)
            if s.set.quasi_static
                solver = FBDF(nlsolve=OrdinaryDiffEqNonlinearSolve.NLNewton(relax=0.4, max_iter=1000))
            else
                solver = FBDF()
            end
            s.sys = s.prob.f.sys
            s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false, save_everystep=false)
            sym_vec = get_unknowns(s)
            s.unknowns_vec = zeros(SimFloat, length(sym_vec))
            generate_getters!(s, sym_vec)
        end
        prn && @info "Initialized integrator in $t seconds"
    end

    init_unknowns_vec!(s, s.point_system, s.unknowns_vec, init_Q_b_w, init_kite_pos)
    s.set_unknowns(s.integrator, s.unknowns_vec)
    set_t!(s.integrator, 0.0)
    OrdinaryDiffEqCore.reinit!(s.integrator, s.integrator.u; reinit_dae=true)
    linearize_vsm!(s)
    return nothing
end

function generate_getters!(s, sym_vec)
    sys = s.sys
    c = collect

    set_set_values = setp(sys, sys.set_values)
    set_measure = setp(sys, sys.measured_wind_dir_gnd)
    set_vsm = setp(sys, c.([
        sys.last_x,
        sys.last_y,
        sys.vsm_jac,
    ]))
    set_unknowns = setu(sys, sym_vec)

    get_state = getu(sys, 
        [c(sys.pos), c(sys.acc), c(sys.Q_b_w), sys.elevation, sys.azimuth, sys.course, sys.heading_y, 
        c(sys.e_x), c(sys.tether_vel), c(sys.twist_angle), c(sys.kite_vel)]
    )
    get_y = getu(sys, sys.y)

    s.set_set_values = (integ, val) -> set_set_values(integ, val)
    s.set_measure = (integ, val) -> set_measure(integ, val)
    s.set_vsm = (integ, val) -> set_vsm(integ, val)
    s.set_unknowns = (integ, val) -> set_unknowns(integ, val)

    s.get_state = (integ) -> get_state(integ)
    s.get_y = (integ) -> get_y(integ)
    nothing
end

function linearize_vsm!(s::RamAirKite)
    y = s.get_y(s.integrator)
    jac, x = VortexStepMethod.linearize(
        s.vsm_solver, 
        s.aero, 
        y;
        va_idxs=1:3, 
        omega_idxs=4:6,
        theta_idxs=7:6+length(s.point_system.groups),
        moment_frac=s.bridle_fracs[s.point_system.groups[1].fixed_index])
    s.set_vsm(s.integrator, [x, y, jac])
    nothing
end

function next_step!(s::RamAirKite, set_values=nothing; measure::Union{Measurement, Nothing}=nothing, dt=1/s.set.sample_freq, vsm_interval=1)
    if (!isnothing(set_values)) 
        s.set_set_values(s.integrator, set_values)
    end
    if (!isnothing(measure))
        s.set_measure(s.integrator, measure.wind_dir_gnd)
    end
    if vsm_interval != 0 && s.iter % vsm_interval == 0
        linearize_vsm!(s)
    end
    
    s.t_0 = s.integrator.t
    steptime = @elapsed OrdinaryDiffEqCore.step!(s.integrator, dt, true)
    if !successful_retcode(s.integrator.sol)
        println("Return code for solution: ", s.integrator.sol.retcode)
    end
    @assert successful_retcode(s.integrator.sol)
    s.iter += 1
    s.integrator.t, steptime
end

function get_prob_name(set::Settings; precompile=false)
    suffix = ""
    if precompile
        suffix = ".default"
    end
    if set.quasi_static
        return "prob_static_" * string(set.segments) * "_seg.bin" * suffix
    else
        return "prob_dynamic_" * string(set.segments) * "_seg.bin" * suffix
    end
end

"""
Calculate and return the angle of attack in rad
"""
function calc_aoa(s::RamAirKite)
    alpha_array = s.vsm_solver.sol.alpha_array
    middle = length(alpha_array) ÷ 2
    if iseven(length(alpha_array))
        return 0.5alpha_array[middle] + 0.5alpha_array[middle+1]
    else
        return alpha_array[middle+1]
    end
end

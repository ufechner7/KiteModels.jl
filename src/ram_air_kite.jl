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

"""
    mutable struct RamAirKite{S, V, P} <: AbstractKiteModel

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
    "Linearization problem of the mtk model"
    lin_prob::Union{ModelingToolkit.LinearizationProblem, Nothing} = nothing
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
    defaults::Vector{Pair{Num, Real}} = Pair{Num, Real}[]
    guesses::Vector{Pair{Num, Real}} = Pair{Num, Real}[]

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

function RamAirKite(set::Settings)
    wing = RamAirWing(set; prn=false)
    aero = BodyAerodynamics([wing])
    vsm_solver = Solver(aero; solver_type=NONLIN, atol=2e-8, rtol=2e-8)
    point_system = PointMassSystem(set, wing)
    return RamAirKite(set, aero, vsm_solver, point_system)
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
    init_sim!(s::RamAirKite, measure::Measurement; prn=true, precompile=false) -> Nothing

Initialize a kite power system model. 

If a serialized model exists for the current configuration, it will load that model
and only update the state variables. Otherwise, it will create a new model from scratch.

# Fast path (serialized model exists):
1. Loads existing ODEProblem from disk
2. Calls `reinit!` to update state variables
3. Sets up integrator with current measurements

# Slow path (no serialized model):
1. Creates symbolic MTK system with all equations
2. Simplifies system equations
3. Creates ODEProblem and serializes to disk
4. Proceeds with fast path

# Arguments
- `s::RamAirKite`: The kite system state object  
- `measure::Measurement`: Initial state measurements
- `prn::Bool=true`: Whether to print progress information
- `precompile::Bool=false`: Whether to build problem for precompilation

# Returns
`Nothing`
"""
function init_sim!(s::RamAirKite, measure::Measurement;
    solver=ifelse(s.set.quasi_static, FBDF(nlsolve=OrdinaryDiffEqNonlinearSolve.NLNewton(relax=0.4, max_iter=1000)), FBDF()), 
    adaptive=true, prn=true, precompile=false, lin_sys=false, remake=false, reload=false
)
    function init(s, measure)
        init_Q_b_w, R_b_w = measure_to_q(measure)
        init_kite_pos = init!(s.point_system, s.set, R_b_w)

        init_va_b = R_b_w' * [s.set.v_wind, 0., 0.]
        
        sys, defaults, guesses, inputs = create_sys!(s, s.point_system, measure; init_Q_b_w, init_kite_pos, init_va_b, lin_sys)
        prn && @info "Simplifying the system"
        prn ? (@time sys = structural_simplify(sys; additional_passes=[ModelingToolkit.IfLifting])) :
            (sys = structural_simplify(sys; additional_passes=[ModelingToolkit.IfLifting]))
        s.sys = sys
        dt = SimFloat(1/s.set.sample_freq)
        if prn
            @info "Creating ODEProblem"
            @time s.prob = ODEProblem(s.sys, defaults, (0.0, dt); guesses, initializealg=CheckInit())
        else
            s.prob = ODEProblem(s.sys, defaults, (0.0, dt); guesses, initializealg=CheckInit())
        end
        serialize(prob_path, s.prob)
        s.integrator = nothing
    end
    prob_path = joinpath(KiteUtils.get_data_path(), get_prob_name(s.set; precompile))
    if !ispath(prob_path) || remake
        init(s, measure)
    end
    success = reinit!(s, measure; solver, adaptive, precompile, reload)
    if !success
        rm(prob_path)
        @info "Rebuilding the system. This can take some minutes..."
        init(s, measure)
        reinit!(s, measure; precompile, prn)
    end
    return nothing
end

function init_sim!(::RamAirKite; prn=true)
    throw(ArgumentError("Use the function init_sim!(s::RamAirKite, measure::Measurement) instead."))
end

"""
    reinit!(s::RamAirKite; prn=true, precompile=false) -> Nothing

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
function reinit!(
    s::RamAirKite, measure::Measurement; 
    solver=ifelse(s.set.quasi_static, FBDF(nlsolve=OrdinaryDiffEqNonlinearSolve.NLNewton(relax=0.4, max_iter=1000)), FBDF()),
    adaptive=true,
    prn=true, 
    reload=true, 
    precompile=false
)
    isnothing(s.point_system) && throw(ArgumentError("PointMassSystem not defined"))

    init_Q_b_w, R_b_w = measure_to_q(measure)
    init_kite_pos = init!(s.point_system, s.set, R_b_w)
    
    if isnothing(s.prob) || reload
        prob_path = joinpath(KiteUtils.get_data_path(), get_prob_name(s.set; precompile))
        !ispath(prob_path) && throw(ArgumentError("$prob_path not found. Run init_sim!(s::RamAirKite) first."))
        try
            s.prob = deserialize(prob_path)
        catch e
            @warn "Failure to deserialize $prob_path !"
            return false
        end
    end
    if isnothing(s.integrator) || !successful_retcode(s.integrator.sol) || reload
        t = @elapsed begin
            dt = SimFloat(1/s.set.sample_freq)
            s.sys = s.prob.f.sys
            s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; 
                adaptive, dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false, save_everystep=false)
            sym_vec = get_unknowns(s)
            s.unknowns_vec = zeros(SimFloat, length(sym_vec))
            generate_getters!(s, sym_vec)
        end
        prn && @info "Initialized integrator in $t seconds"
    end

    init_unknowns_vec!(s, s.point_system, s.unknowns_vec, init_Q_b_w, init_kite_pos)
    s.set_unknowns(s.integrator, s.unknowns_vec)
    OrdinaryDiffEqCore.reinit!(s.integrator, s.integrator.u; reinit_dae=true)
    linearize_vsm!(s)
    return true
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
        [c(sys.pos), c(sys.acc), c(sys.Q_b_w), sys.elevation, sys.azimuth, sys.course, sys.heading_x, 
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
        moment_frac=0.0)
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
    ver = "$(VERSION.major).$(VERSION.minor)"
    if precompile
        suffix = ".default"
    end
    dynamics_type = ifelse(set.quasi_static, "static", "dynamic")
    return "prob_$(ver)_$(set.physical_model)_$(dynamics_type)_$(set.segments)_seg.bin$suffix"
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

function init_unknowns_vec!(
    s::RamAirKite, 
    system::PointMassSystem, 
    vec::Vector{SimFloat},
    init_Q_b_w,
    init_kite_pos;
    non_observed=true
)
    !s.set.quasi_static && non_observed && (length(vec) != length(s.integrator.u)) && 
        throw(ArgumentError("Unknowns of length $(length(s.integrator.u)) but vector provided of length $(length(vec))"))
        
    @unpack points, groups, segments, pulleys, winches = system
    vec_idx = 1
    
    if non_observed
        for point in points
            if point.type == DYNAMIC
                for i in 1:3
                    vec[vec_idx] = point.pos_w[i]
                    vec_idx += 1
                end
                for i in 1:3 # TODO: add speed to init
                    vec[vec_idx] = 0.0
                    vec_idx += 1
                end
            end
        end
        
        for pulley in pulleys
            if pulley.type == DYNAMIC
                vec[vec_idx] = segments[pulley.segments[1]].l0
                vec_idx += 1
                vec[vec_idx] = 0
                vec_idx += 1
            end
        end
    end

    for group in groups
        if group.type == DYNAMIC
            vec[vec_idx] = 0
            vec_idx += 1
            vec[vec_idx] = 0
            vec_idx += 1
        end
    end

    for winch in winches
        vec[vec_idx] = winch.tether_length
        vec_idx += 1
        vec[vec_idx] = 0
        vec_idx += 1
    end

    for i in 1:4
        vec[vec_idx] = init_Q_b_w[i]
        vec_idx += 1
    end
    for i in 1:3
        vec[vec_idx] = 0
        vec_idx += 1
    end
    for i in 1:3
        vec[vec_idx] = init_kite_pos[i]
        vec_idx += 1
    end
    for i in 1:3
        vec[vec_idx] = 0
        vec_idx += 1
    end
    non_observed && (vec_idx-1 != length(vec)) && 
        throw(ArgumentError("Unknowns vec is of length $(length(vec)) but the last index is $(vec_idx-1)"))
    nothing
end

function get_unknowns(s::RamAirKite; simple=false)
    vec = Num[]
    vec = get_stiff_unknowns(s, vec)
    vec = get_nonstiff_unknowns(s, vec; simple)
    !s.set.quasi_static && (length(vec) != length(s.integrator.u)) &&
        throw(ArgumentError("Integrator unknowns of length $(length(s.integrator.u)) should equal vec of length $(length(vec))"))
    return vec
end

function get_stiff_unknowns(s, vec=Num[])
    @unpack points, groups, segments, pulleys, winches = s.point_system
    for point in points
        for i in 1:3
            point.type == DYNAMIC && push!(vec, s.sys.pos[i, point.idx])
        end
        for i in 1:3 # TODO: add speed to init
            point.type == DYNAMIC && push!(vec, s.sys.vel[i, point.idx])
        end
    end
    for pulley in pulleys
        pulley.type == DYNAMIC && push!(vec, s.sys.pulley_l0[pulley.idx])
        pulley.type == DYNAMIC && push!(vec, s.sys.pulley_vel[pulley.idx])
    end
    return vec
end

function get_nonstiff_unknowns(s, vec=Num[]; simple=false)
    @unpack points, groups, segments, pulleys, winches = s.point_system
    if !simple
        for group in groups
            group.type == DYNAMIC && push!(vec, s.sys.free_twist_angle[group.idx])
            group.type == DYNAMIC && push!(vec, s.sys.twist_ω[group.idx])
        end
    else
        for i in 1:2
            push!(vec, s.sys.simple_twist_angle[i])
            push!(vec, s.sys.simple_twist_ω[i])
        end
    end
    for winch in winches
        push!(vec, s.sys.tether_length[winch.idx])
        push!(vec, s.sys.tether_vel[winch.idx])
    end
    [push!(vec, s.sys.Q_b_w[i]) for i in 1:4]
    [push!(vec, s.sys.ω_b[i]) for i in 1:3]
    [push!(vec, s.sys.kite_pos[i]) for i in 1:3]
    [push!(vec, s.sys.kite_vel[i]) for i in 1:3]
    return vec
end


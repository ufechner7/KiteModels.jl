# Copyright (c) 2024, 2025 Bart van de Lint and Uwe Fechner
# SPDX-License-Identifier: MIT

@with_kw mutable struct SerializedModel
    set_hash::Vector{UInt8}
    "Reference to the geometric wing model"
    vsm_wings::Vector{VortexStepMethod.RamAirWing}
    "Reference to the aerodynamic wing model"
    vsm_aeros::Vector{VortexStepMethod.BodyAerodynamics}
    "Reference to the VSM aerodynamics solver"
    vsm_solvers::Vector{VortexStepMethod.Solver}
    sys_struct_hash::Vector{UInt8}
    "Reference to the atmospheric model as implemented in the package AtmosphericModels"
    am::AtmosphericModel = AtmosphericModel()
    "Simplified system of the mtk model"
    sys::Union{ModelingToolkit.ODESystem, Nothing} = nothing
    "Unsimplified system of the mtk model"
    full_sys::Union{ModelingToolkit.ODESystem, Nothing} = nothing
    "Linearization function of the mtk model"
    lin_prob::Union{ModelingToolkit.LinearizationProblem, Nothing} = nothing
    "ODE function of the mtk model"
    prob::Union{OrdinaryDiffEqCore.ODEProblem, Nothing} = nothing

    unknowns_vec::Vector{SimFloat} = zeros(SimFloat, 3)
    defaults::Vector{Pair} = Pair[]
    guesses::Vector{Pair} = Pair[]

    set_psys::Union{Function, Nothing}             = nothing
    set_set_values::Union{Function, Nothing}       = nothing
    set_set::Union{Function, Nothing}             = nothing
    set_vsm::Union{Function, Nothing}              = nothing
    set_unknowns::Union{Function, Nothing}         = nothing
    set_nonstiff::Union{Function, Nothing}         = nothing
    set_lin_vsm::Union{Function, Nothing}          = nothing
    set_lin_set_values::Union{Function, Nothing}   = nothing
    set_lin_unknowns::Union{Function, Nothing}     = nothing
    set_stabilize::Union{Function, Nothing}        = nothing
    
    get_vsm::Union{Function, Nothing}              = nothing
    get_set_values::Union{Function, Nothing}       = nothing
    get_unknowns::Union{Function, Nothing}         = nothing
    get_wing_state::Union{Function, Nothing}       = nothing
    get_winch_state::Union{Function, Nothing}      = nothing
    get_point_state::Union{Function, Nothing}      = nothing
    get_y::Union{Function, Nothing}                = nothing
    get_unstretched_length::Union{Function, Nothing} = nothing
    get_tether_length::Union{Function, Nothing}    = nothing
    get_wing_pos::Union{Function, Nothing}         = nothing
    get_winch_force::Union{Function, Nothing}      = nothing
    get_spring_force::Union{Function, Nothing}     = nothing
    get_stabilize::Union{Function, Nothing}        = nothing
    get_pos::Union{Function, Nothing}              = nothing
end

"""
    mutable struct SymbolicAWEModel{S, V, P} <: AbstractKiteModel

State of the kite power system, using a quaternion kite model and three steering lines to the ground. Parameters:
- S: Scalar type, e.g. SimFloat
  In the documentation mentioned as Any, but when used in this module it is always SimFloat and not Any.
- V: Vector type, e.g. KVec3
- P: number of tether points of the system, 3segments+3
Normally a user of this package will not have to access any of the members of this type directly,
use the input and output functions instead.

$(TYPEDFIELDS)
"""
@with_kw mutable struct SymbolicAWEModel <: AbstractKiteModel
    "Reference to the settings struct"
    set::Settings
    "Reference to the point mass system with points, segments, pulleys and tethers"
    sys_struct::SystemStructure
    serialized_model::SerializedModel
    integrator::Union{OrdinaryDiffEqCore.ODEIntegrator, Nothing} = nothing
    "relative start time of the current time interval"
    t_0::SimFloat = 0.0
    "Number of solve! calls"
    iter::Int64 = 0
    t_vsm::SimFloat  = zero(SimFloat)
    t_step::SimFloat = zero(SimFloat)
    set_tether_length::Vector{SimFloat} = zeros(SimFloat, 3)
end

function Base.getproperty(sam::SymbolicAWEModel, sym::Symbol)
    if hasfield(SymbolicAWEModel, sym)
        getfield(sam, sym)
    else
        getproperty(getfield(sam, :serialized_model), sym)
    end
end
function Base.setproperty!(sam::SymbolicAWEModel, sym::Symbol, val)
    if hasfield(SymbolicAWEModel, sym)
        setfield!(sam, sym, val)
    else
        serialized_model = getfield(sam, :serialized_model)
        setproperty!(serialized_model, sym, val)
    end
end

"""
    SymbolicAWEModel(set::Settings, sys_struct::SystemStructure, 
                     vsm_aeros::Vector{<:BodyAerodynamics}=BodyAerodynamics[], 
                     vsm_solvers::Vector{<:VortexStepMethod.Solver}=VortexStepMethod.Solver[])

Constructs a SymbolicAWEModel that can generate ModelingToolkit equations
from the discrete mass-spring-damper representation defined in the [`SystemStructure`](@ref).
The aerodynamic models provide forces and moments acting on wing components.

# Arguments
- `set::Settings`: Configuration parameters (see [KiteUtils.Settings](https://ufechner7.github.io/KiteUtils.jl/stable/types/#KiteUtils.Settings))
- `sys_struct::SystemStructure`: Physical system definition with points, segments, groups, etc.
- `vsm_aeros::Vector{<:BodyAerodynamics}=BodyAerodynamics[]`: Aerodynamic models for each wing
- `vsm_solvers::Vector{<:VortexStepMethod.Solver}=VortexStepMethod.Solver[]`: VSM solvers for aerodynamic calculations

# Returns
- `SymbolicAWEModel`: Model ready for symbolic equation generation via [`init_sim!`](@ref)

# Example
```julia
# Create wing geometry and aerodynamics
set = se()
wing = RamAirWing(set)
aero = BodyAerodynamics([wing])
solver = Solver(aero; solver_type=NONLIN)

# Create system structure
sys_struct = SystemStructure(set, wing)

# Create symbolic model
model = SymbolicAWEModel(set, sys_struct, [aero], [solver])
```
"""
function SymbolicAWEModel(
    set::Settings, 
    sys_struct::SystemStructure,
    vsm_aeros::Vector{<:BodyAerodynamics}=BodyAerodynamics[], 
    vsm_solvers::Vector{<:VortexStepMethod.Solver}=VortexStepMethod.Solver[]
)
    vsm_wings = [aero.wings[1] for aero in vsm_aeros]
    set_hash = get_set_hash(set)
    sys_struct_hash = get_sys_struct_hash(sys_struct)
    serialized_model = SerializedModel(; set_hash, sys_struct_hash, vsm_wings, vsm_aeros, vsm_solvers)
    return SymbolicAWEModel(; set, sys_struct, serialized_model)
end

"""
    SymbolicAWEModel(set::Settings)

Constructs a default SymbolicAWEModel with automatically generated components.

This convenience constructor creates a complete AWE model using default configurations:
- Generates a ram-air wing from settings
- Creates aerodynamic model and VSM solver
- Builds system structure based on the wing geometry
- Assembles everything into a ready-to-use symbolic model

# Arguments
- `set::Settings`: Configuration parameters (see [KiteUtils.Settings](https://ufechner7.github.io/KiteUtils.jl/stable/types/#KiteUtils.Settings))

# Returns
- `SymbolicAWEModel`: Model ready for symbolic equation generation via [`init_sim!`](@ref)

# Example
```julia
set = se()  # Load default settings
model = SymbolicAWEModel(set)

# Initialize and run simulation
init_sim!(model)
for i in 1:1000
    next_step!(model)
end
```
"""
function SymbolicAWEModel(set::Settings)
    wing = RamAirWing(set; prn=false)
    aero = BodyAerodynamics([wing])
    vsm_solver = Solver(aero; solver_type=NONLIN, atol=2e-8, rtol=2e-8)
    sys_struct = SystemStructure(set, wing)
    return SymbolicAWEModel(set, sys_struct, [aero], [vsm_solver])
end

function update_sys_state!(ss::SysState, s::SymbolicAWEModel, zoom=1.0)
    isnothing(s.integrator) && error("run init_sim!(s) first")
    ss.time = s.integrator.t # Use integrator time

    # Get the state vectors from the integrator
    if !isnothing(s.get_winch_state)
        nw = length(s.sys_struct.winches)
        set_values, tether_length, tether_vel, winch_force = s.get_winch_state(s.integrator)
        ss.l_tether[1:nw] .= tether_length
        ss.v_reelout[1:nw] .= tether_vel
        ss.force[1:nw] .= winch_force
        ss.set_torque[1:nw] .= set_values
    end
    if !isnothing(s.get_wing_state)
        Q_b_w, elevation, azimuth, course, heading, twist_angle, wing_vel, aero_force_b, 
            aero_moment_b, turn_rate, va_wing_b, wind_vel_wing = s.get_wing_state(s.integrator)
        ss.orient .= Q_b_w[1, :]
        ss.turn_rates .= turn_rate[1, :]
        ss.elevation = elevation[1]
        ss.azimuth = azimuth[1]
        # Depower and Steering from twist angles
        num_groups = length(s.sys_struct.wings[1].group_idxs)
        ss.twist_angles[1:num_groups] .= twist_angle[1:num_groups]
        ss.depower = rad2deg(mean(ss.twist_angles)) # Average twist for depower
        ss.steering = rad2deg(ss.twist_angles[num_groups] - ss.twist_angles[1])
        ss.heading = heading[1] # Use heading from MTK model
        ss.course = course[1]
        # Apparent Wind and Aerodynamics
        ss.v_app = norm(va_wing_b[1, :])
        ss.v_wind_kite .= wind_vel_wing[1, :]
        # Calculate AoA and Side Slip from apparent wind in body frame
        # AoA: angle between v_app projected onto xz-plane and x-axis
        # Side Slip: angle between v_app and the xz-plane
        if ss.v_app > 1e-6 # Avoid division by zero
            ss.AoA = atan(va_wing_b[1, 3], va_wing_b[1, 1])
            ss.side_slip = asin(va_wing_b[1, 2] / ss.v_app)
        else
            ss.AoA = 0.0
            ss.side_slip = 0.0
        end
        ss.aero_force_b .= aero_force_b[1, :]
        ss.aero_moment_b .= aero_moment_b[1, :]
        ss.vel_kite .= wing_vel[1, :]
        # Calculate Roll, Pitch, Yaw from Quaternion
        q = Q_b_w[1, :]
        # roll (x-axis rotation)
        sinr_cosp = 2 * (q[1] * q[2] + q[3] * q[4])
        cosr_cosp = 1 - 2 * (q[2] * q[2] + q[3] * q[3])
        ss.roll = atan(sinr_cosp, cosr_cosp)
        # pitch (y-axis rotation)
        sinp = 2 * (q[1] * q[3] - q[4] * q[2])
        if abs(sinp) >= 1
            ss.pitch = copysign(pi / 2, sinp) # use 90 degrees if out of range
        else
            ss.pitch = asin(sinp)
        end
        # yaw (z-axis rotation)
        siny_cosp = 2 * (q[1] * q[4] + q[2] * q[3])
        cosy_cosp = 1 - 2 * (q[3] * q[3] + q[4] * q[4])
        ss.yaw = atan(siny_cosp, cosy_cosp)
    end

    pos, acc, wind_vec_gnd = s.get_point_state(s.integrator)
    P = length(s.sys_struct.points)
    for i in 1:P
        ss.X[i] = pos[1, i] * zoom
        ss.Y[i] = pos[2, i] * zoom
        ss.Z[i] = pos[3, i] * zoom
    end
    ss.acc = norm(acc) # Use the norm of the wing's acceleration vector
    ss.v_wind_gnd .= wind_vec_gnd
    nothing
end

function SysState(s::SymbolicAWEModel, zoom=1.0)
    ss = SysState{length(s.sys_struct.points)}()
    update_sys_state!(ss, s, zoom)
    ss
end

"""
    init_sim!(s::SymbolicAWEModel; solver=nothing, adaptive=true, prn=true, 
              precompile=false, remake=false, reload=false, 
              lin_outputs=Num[]) -> OrdinaryDiffEqCore.ODEIntegrator

Initialize a kite power system model. 

If a serialized model exists for the current configuration, it will load that model
and only update the state variables. Otherwise, it will create a new model from scratch.

# Fast path (serialized model exists):
1. Loads existing ODEProblem from disk
2. Calls `reinit!` to update state variables
3. Sets up integrator with initial settings

# Slow path (no serialized model):
1. Creates symbolic MTK system with all equations
2. Simplifies system equations
3. Creates ODEProblem and serializes to disk
4. Proceeds with fast path

# Arguments
- `s::SymbolicAWEModel`: The kite system state object  

# Keyword arguments
- `solver`: Solver algorithm to use. If `nothing`, defaults to `FBDF()` or `QNDF()` based on `s.set.solver`.
- `adaptive::Bool=true`: Whether to use adaptive time stepping.
- `prn::Bool=true`: Whether to print progress information.
- `precompile::Bool=false`: Whether to build problem for precompilation.
- `remake::Bool=false`: If true, forces the system to be rebuilt, even if a serialized model exists.
- `reload::Bool=false`: If true, forces the system to reload the serialized model from disk.
- `lin_outputs::Vector{Num}=Num[]`: List of symbolic variables for which to generate a linearization function.

# Returns
- `integrator::OrdinaryDiffEqCore.ODEIntegrator`: The initialized ODE integrator.
"""
function init_sim!(s::SymbolicAWEModel; 
    solver=nothing, adaptive=true, prn=true, 
    precompile=false, remake=false, reload=false, 
    lin_outputs=Num[]
)
    if isnothing(solver)
        solver = if s.set.solver == "FBDF"
            if s.set.quasi_static
                FBDF(nlsolve=OrdinaryDiffEqNonlinearSolve.NLNewton(relax=s.set.relaxation))
            else
                FBDF()
            end
        elseif s.set.solver == "QNDF"
            @warn "This solver is not tested."
            QNDF()
        else
            error("Unavailable solver for SymbolicAWEModel: $(s.set.solver).")
        end
    end
    function init(s)
        init_Q_b_w, R_b_w, init_va_b = initial_orient(s)
        init!(s.sys_struct, s.set)
        
        inputs = create_sys!(s, s.sys_struct; init_va_b)
        prn && @info "Simplifying the system"
        prn ? (@time (sys, _) = structural_simplify(s.full_sys, (inputs, []))) :
            ((sys, _) = structural_simplify(s.full_sys, (inputs, [])))
        s.sys = sys
        dt = SimFloat(1/s.set.sample_freq)
        if prn
            @info "Creating ODEProblem"
            @time s.prob = ODEProblem(s.sys, s.defaults, (0.0, dt); s.guesses)
        else
            s.prob = ODEProblem(s.sys, s.defaults, (0.0, dt); s.guesses)
        end
        if length(lin_outputs) > 0
            lin_fun, _ = linearization_function(s.full_sys, [inputs...], lin_outputs; op=s.defaults, guesses=s.guesses)
            s.lin_prob = LinearizationProblem(lin_fun, 0.0)
        end
        sym_vec = get_unknowns(s.sys_struct, s.sys)
        s.unknowns_vec = zeros(SimFloat, length(sym_vec))
        generate_getters!(s, sym_vec)
        s.set_hash = get_set_hash(s.set)
        s.sys_struct_hash = get_sys_struct_hash(s.sys_struct)
        serialize(model_path, s.serialized_model)
        s.integrator = nothing
        return nothing
    end
    model_path = joinpath(KiteUtils.get_data_path(), get_model_name(s.set; precompile))
    if !ispath(model_path) || remake
        init(s)
    end
    _, success = reinit!(s, solver; adaptive, precompile, reload, lin_outputs, prn)
    if !success
        rm(model_path)
        @info "Rebuilding the system. This can take some minutes..."
        init(s)
        reinit!(s, solver; adaptive, precompile, lin_outputs, prn, reload=true)
    end
    return s.integrator
end

function linearize(s::SymbolicAWEModel; set_values=s.get_set_values(s.integrator))
    isnothing(s.lin_prob) && error("Run init_sim! with remake=true and lin_outputs=...")
    s.set_lin_vsm(s.lin_prob, s.get_vsm(s.integrator))
    s.set_lin_set_values(s.lin_prob, set_values)
    s.set_lin_unknowns(s.lin_prob, s.get_unknowns(s.integrator))
    return solve(s.lin_prob)
end

"""
    reinit!(s::SymbolicAWEModel, solver; prn=true, precompile=false) -> Nothing

Reinitialize an existing kite power system model with new state values.
The new state is coming from the init section of the settings, stored
in the struct `s.set`.

This function performs the following operations:
1. If no integrator exists yet:
   - Loads a serialized ODEProblem from disk
   - Initializes a new ODE integrator 
   - Generates getter/setter functions for the system
2. Converts initial settings to quaternion orientation
3. Initializes the point mass system with new positions
4. Sets initial values for all state variables
5. Reinitializes the ODE integrator
6. Updates the linearized aerodynamic model

This is more efficient than `init!` as it reuses the existing model structure
and only updates the state variables to match the current initial settings.

# Arguments
- `s::SymbolicAWEModel`: The kite power system state object
- `solver`: The solver to be used
- `prn::Bool=true`: Whether to print progress information

# Returns
- `Nothing`

# Throws
- `ArgumentError`: If no serialized problem exists (run `init_sim!` first)
"""
function reinit!(
    s::SymbolicAWEModel,
    solver;
    adaptive=true,
    prn=true, 
    reload=true, 
    precompile=false,
    lin_outputs=Num[]
)
    isnothing(s.sys_struct) && error("SystemStructure not defined")

    # init_Q_b_w, R_b_w, init_va_b = initial_orient(s)
    
    if isnothing(s.prob) || reload
        model_path = joinpath(KiteUtils.get_data_path(), get_model_name(s.set; precompile))
        !ispath(model_path) && error("$model_path not found. Run init_sim!(s::SymbolicAWEModel) first.")
        try
            s.serialized_model = deserialize(model_path)
        catch e
            @warn "Failure to deserialize $model_path !"
            return s.integrator, false
        end
        if length(lin_outputs) > 0 && isnothing(s.lin_prob) 
            @warn "lin_prob is nothing."
            return s.integrator, false
        elseif (get_set_hash(s.set) != s.serialized_model.set_hash)
            @warn "The Settings have changed."
            return s.integrator, false
        elseif (get_sys_struct_hash(s.sys_struct) != s.serialized_model.sys_struct_hash)
            @warn "The SystemStructure has changed."
            return s.integrator, false
        end
    end
    if isnothing(s.integrator) || !successful_retcode(s.integrator.sol) || reload
        t = @elapsed begin
            dt = SimFloat(1/s.set.sample_freq)
            s.sys = s.prob.f.sys
            s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; 
                adaptive, dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false, save_everystep=false)
            !s.set.quasi_static && (length(s.unknowns_vec) != length(s.integrator.u)) &&
                error("sam.integrator unknowns of length $(length(s.integrator.u)) should equal sam.unknowns_vec of length $(length(s.unknowns_vec)).
                    Maybe you forgot to run init_sim!(model; remake=true)?")
        end
        prn && @info "Initialized integrator in $t seconds"
    end

    init!(s.sys_struct, s.set)
    init_unknowns_vec!(s, s.sys_struct, s.unknowns_vec)
    s.set_unknowns(s.integrator, s.unknowns_vec)
    s.set_psys(s.integrator, s.sys_struct)
    OrdinaryDiffEqCore.reinit!(s.integrator, s.integrator.u; reinit_dae=true)
    linearize_vsm!(s)
    return s.integrator, true
end

function generate_getters!(s, sym_vec)
    sys = s.sys
    c = collect
    @unpack wings, winches, tethers = s.sys_struct

    if length(wings) > 0
        vsm_sym = c.([
            sys.last_x,
            sys.last_y,
            sys.vsm_jac,
        ])
        get_vsm = getp(sys, vsm_sym)
        s.get_vsm = (integ) -> get_vsm(integ)
        get_y = getu(sys, sys.y)
        s.get_y = (integ) -> get_y(integ)
        get_wing_state = getu(sys, c.([
            sys.Q_b_w,           # Orientation quaternion
            sys.elevation,       # Elevation angle
            sys.azimuth,         # Azimuth angle
            sys.course,          # Course angle
            sys.heading,         # Heading angle (based on body x-axis projection)
            sys.twist_angle,     # Twist angle per group
            sys.wing_vel,        # Kite center velocity vector (world frame)
            sys.aero_force_b,    # Aerodynamic force (body frame)
            sys.aero_moment_b,   # Aerodynamic moment (body frame)
            sys.turn_rate,       # Angular velocity (body frame)
            sys.va_wing_b,       # Apparent wind velocity (body frame)
            sys.wind_vel_wing    # Wind vector at wing height (world frame)
        ]))
        s.get_wing_state = (integ) -> get_wing_state(integ)
        get_wing_pos = getu(sys, sys.wing_pos)
        s.get_wing_pos = (integ) -> get_wing_pos(integ)

        set_vsm = setp(sys, vsm_sym)
        s.set_vsm = (integ, val) -> set_vsm(integ, val)
        if !isnothing(s.lin_prob)
            set_lin_vsm = setp(s.lin_prob, vsm_sym)
            s.set_lin_vsm = (lin_prob, val) -> set_lin_vsm(lin_prob, val)
        end
    end

    if length(winches) > 0
        get_winch_state = getu(sys, c.([
             sys.set_values,
             sys.tether_length,   # Unstretched length per winch
             sys.tether_vel,      # Reeling velocity per winch
             sys.winch_force,     # Force at winch connection point per winch
        ]))
        s.get_winch_state = (integ) -> get_winch_state(integ)
        get_set_values = getp(sys, sys.set_values)
        s.get_set_values = (integ) -> get_set_values(integ)
        get_tether_length = getu(sys, sys.tether_length)
        s.get_tether_length = (integ) -> get_tether_length(integ)
        get_winch_force = getu(sys, sys.winch_force)
        s.get_winch_force = (integ) -> get_winch_force(integ)

        set_set_values = setp(sys, sys.set_values)
        s.set_set_values = (integ, val) -> set_set_values(integ, val)
        if !isnothing(s.lin_prob)
            set_lin_set_values = setp(s.lin_prob, sys.set_values)
            s.set_lin_set_values = (lin_prob, val) -> set_lin_set_values(lin_prob, val)
        end
    end

    if length(tethers) > 0
        get_unstretched_length = getu(sys, sys.tether_length)
        s.get_unstretched_length = (integ) -> get_unstretched_length(integ)
    end

    if length(winches) + length(wings) > 0
        set_stabilize = setp(sys, sys.stabilize)
        s.set_stabilize = (integ, val) -> set_stabilize(integ, val)

        get_stabilize = getp(sys, sys.stabilize)
        s.get_stabilize = (integ) -> get_stabilize(integ)
    end
    
    set_psys = setp(sys, sys.psys)
    s.set_psys = (integ, val) -> set_psys(integ, val)
    set_set = setp(sys, sys.pset)
    s.set_set = (integ, val) -> set_set(integ, val)
    set_unknowns = setu(sys, sym_vec)
    s.set_unknowns = (integ, val) -> set_unknowns(integ, val)
    set_nonstiff = setu(sys, get_nonstiff_unknowns(s.sys_struct, s.sys))
    s.set_nonstiff = (integ, val) -> set_nonstiff(integ, val)
    
    get_unknowns = getu(sys, sym_vec)
    s.get_unknowns = (integ) -> get_unknowns(integ)
    get_point_state = getu(sys, c.([
         sys.pos,             # Particle positions
         sys.acc,             # Kite center acceleration vector (world frame)
         sys.wind_vec_gnd,    # Ground wind vector (world frame)
    ]))
    s.get_point_state = (integ) -> get_point_state(integ)
    get_spring_force = getu(sys, sys.spring_force)
    s.get_spring_force = (integ) -> get_spring_force(integ)
    get_pos = getu(sys, sys.pos)
    s.get_pos = (integ) -> get_pos(integ)
    
    if !isnothing(s.lin_prob)
        set_lin_unknowns = setu(s.lin_prob, Initial.(sym_vec))
        s.set_lin_unknowns = (lin_prob, val) -> set_lin_unknowns(lin_prob, val)
    end
    nothing
end

function linearize_vsm!(s::SymbolicAWEModel, integ=s.integrator)
    @unpack wings, y, x, jac = s.sys_struct
    if length(wings) > 0
        y .= s.get_y(integ)
        for wing in wings
            res = VortexStepMethod.linearize(
                s.vsm_solvers[wing.idx], 
                s.vsm_aeros[wing.idx], 
                y[wing.idx, :];
                va_idxs=1:3, 
                omega_idxs=4:6,
                theta_idxs=7:6+length(s.sys_struct.groups),
                moment_frac=s.sys_struct.groups[1].moment_frac
            )
            jac[wing.idx, :, :] .= res[1]
            x[wing.idx, :] .= res[2]
        end
        s.set_vsm(integ, [x, y, jac])
    end
    nothing
end

"""
    next_step!(s::SymbolicAWEModel, integrator::ODEIntegrator; set_values=nothing, upwind_dir=nothing, dt=1/s.set.sample_freq, vsm_interval=1)

Take a simulation step, using the internal integrator.

This function performs the following steps:
1. Optionally update the set values (control inputs)
2. Optionally update the upwind direction
3. Optionally linearize the VSM (Vortex Step Method) model
4. Step the ODE integrator forward by `dt` seconds
5. Check for a successful return code from the integrator
6. Increment the iteration counter

# Arguments
- `s::SymbolicAWEModel`: The kite power system state object
- `integrator::ODEIntegrator`: The ODE integrator to use

# Keyword Arguments
- `set_values=nothing`: New values for the set variables (control inputs). If `nothing`, the current values are used.
- `dt=1/s.set.sample_freq`: Time step size in seconds. Defaults to the inverse of the sample frequency.
- `vsm_interval=1`: Interval (in number of steps) at which to linearize the VSM model. If 0, the VSM model is not linearized.

# Returns
- `Nothing`
"""
function next_step!(s::SymbolicAWEModel, integrator::OrdinaryDiffEqCore.ODEIntegrator; set_values=nothing, dt=1/s.set.sample_freq, vsm_interval=1)
    !(s.integrator === integrator) && error("The ODEIntegrator doesn't belong to the SymbolicAWEModel")
    next_step!(s; set_values, upwind_dir, dt, vsm_interval)
end

function next_step!(s::SymbolicAWEModel; set_values=nothing, dt=1/s.set.sample_freq, vsm_interval=1)
    if (!isnothing(set_values)) 
        s.set_set_values(s.integrator, set_values)
    end
    if vsm_interval != 0 && s.iter % vsm_interval == 0
        s.t_vsm = @elapsed linearize_vsm!(s)
    end
    
    s.t_0 = s.integrator.t
    s.t_step = @elapsed OrdinaryDiffEqCore.step!(s.integrator, dt, true)
    if !successful_retcode(s.integrator.sol)
        @warn "Return code for solution: $(s.integrator.sol.retcode)"
    end
    @assert successful_retcode(s.integrator.sol)
    s.iter += 1
    return nothing
end

function get_model_name(set::Settings; precompile=false)
    suffix = ""
    ver = "$(VERSION.major).$(VERSION.minor)"
    if precompile
        suffix = ".default"
    end
    dynamics_type = ifelse(set.quasi_static, "static", "dynamic")
    return "model_$(ver)_$(set.physical_model)_$(dynamics_type)_$(set.segments)_seg.bin$suffix"
end

"""
Calculate and return the angle of attack in rad
"""
function calc_aoa(s::SymbolicAWEModel)
    alpha_array = s.vsm_solvers[1].sol.alpha_array
    middle = length(alpha_array) ÷ 2
    if iseven(length(alpha_array))
        return 0.5alpha_array[middle] + 0.5alpha_array[middle+1]
    else
        return alpha_array[middle+1]
    end
end

function init_unknowns_vec!(
    s::SymbolicAWEModel, 
    system::SystemStructure, 
    vec::Vector{SimFloat}
)
    !s.set.quasi_static && (length(vec) != length(s.integrator.u)) && 
        error("Unknowns of length $(length(s.integrator.u)) but vector provided of length $(length(vec))")
        
    @unpack points, groups, segments, pulleys, winches, wings = system
    vec_idx = 1
    
    for point in points
        if point.type == DYNAMIC
            for i in 1:3
                vec[vec_idx] = point.pos_w[i]
                vec_idx += 1
            end
            for i in 1:3
                vec[vec_idx] = point.vel_w[i]
                vec_idx += 1
            end
        end
    end
    for pulley in pulleys
        if pulley.type == DYNAMIC
            vec[vec_idx] = pulley.length
            vec_idx += 1
            vec[vec_idx] = pulley.vel
            vec_idx += 1
        end
    end
    for group in groups
        if group.type == DYNAMIC
            vec[vec_idx] = group.twist
            vec_idx += 1
            vec[vec_idx] = group.twist_vel
            vec_idx += 1
        end
    end
    for winch in winches
        vec[vec_idx] = winch.tether_length
        vec_idx += 1
        vec[vec_idx] = winch.tether_vel
        vec_idx += 1
    end

    for wing in wings
        for i in 1:4
            vec[vec_idx] = wing.orient[i]
            vec_idx += 1
        end
        for i in 1:3
            vec[vec_idx] = wing.angular_vel[i]
            vec_idx += 1
        end
        for i in 1:3
            vec[vec_idx] = wing.pos_w[i]
            vec_idx += 1
        end
        for i in 1:3
            vec[vec_idx] = wing.vel_w[i]
            vec_idx += 1
        end
    end

    (vec_idx-1 != length(vec)) && 
        error("Unknowns vec is of length $(length(vec)) but the last index is $(vec_idx-1)")
    nothing
end

function get_unknowns(sys_struct::SystemStructure, sys::ODESystem)
    vec = Num[]
    @unpack points, groups, segments, pulleys, winches, wings = sys_struct
    for point in points
        for i in 1:3
            point.type == DYNAMIC && push!(vec, sys.pos[i, point.idx])
        end
        for i in 1:3
            point.type == DYNAMIC && push!(vec, sys.vel[i, point.idx])
        end
    end
    for pulley in pulleys
        pulley.type == DYNAMIC && push!(vec, sys.pulley_l0[pulley.idx])
        pulley.type == DYNAMIC && push!(vec, sys.pulley_vel[pulley.idx])
    end
    vec = get_nonstiff_unknowns(sys_struct, sys, vec)
    return vec
end

function get_nonstiff_unknowns(sys_struct::SystemStructure, sys::ODESystem, vec=Num[])
    @unpack points, groups, segments, pulleys, winches, wings = sys_struct
    for group in groups
        group.type == DYNAMIC && push!(vec, sys.free_twist_angle[group.idx])
        group.type == DYNAMIC && push!(vec, sys.twist_ω[group.idx])
    end
    for winch in winches
        push!(vec, sys.tether_length[winch.idx])
        push!(vec, sys.tether_vel[winch.idx])
    end
    for wing in wings
        [push!(vec, sys.Q_b_w[wing.idx, i]) for i in 1:4]
        [push!(vec, sys.ω_b[wing.idx, i]) for i in 1:3]
        [push!(vec, sys.wing_pos[wing.idx, i]) for i in 1:3]
        [push!(vec, sys.wing_vel[wing.idx, i]) for i in 1:3]
    end
    return vec
end

function find_steady_state!(s::SymbolicAWEModel; dt=1/s.set.sample_freq)
    old_state = s.get_stabilize(s.integrator)
    s.set_stabilize(s.integrator, true)
    for _ in 1:1÷dt
        next_step!(s; dt, vsm_interval=1)
    end
    s.set_stabilize(s.integrator, old_state)
    return nothing
end

function initial_orient(s::SymbolicAWEModel)
    set = s.set
    wings = s.sys_struct.wings
    R_b_w = zeros(length(wings), 3, 3)
    Q_b_w = zeros(length(wings), 4)
    init_va_b = zeros(length(wings), 3)
    for wing in wings
        R_cad_body = s.vsm_wings[wing.idx].R_cad_body
        x = [0, 0, -1] # laying flat along x axis
        z = [1, 0, 0] # laying flat along x axis
        x = rotate_around_y(x, -deg2rad(set.elevation))
        z = rotate_around_y(z, -deg2rad(set.elevation))
        x = rotate_around_z(x, deg2rad(set.azimuth))
        z = rotate_around_z(z, deg2rad(set.azimuth))
        R_b_w[wing.idx, :, :] .= R_cad_body' * hcat(x, z × x, z)
        Q_b_w[wing.idx, :] .= rotation_matrix_to_quaternion(R_b_w[wing.idx, :, :])
        init_va_b[wing.idx, :] .= R_b_w[wing.idx, :, :]' * [s.set.v_wind, 0., 0.]
    end
    return Q_b_w, R_b_w, init_va_b
end

unstretched_length(s::SymbolicAWEModel) = s.get_unstretched_length(s.integrator)
tether_length(s::SymbolicAWEModel) = s.get_tether_length(s.integrator)
calc_height(s::SymbolicAWEModel) = s.get_wing_pos(s.integrator)[3]
winch_force(s::SymbolicAWEModel) = s.get_winch_force(s.integrator)
spring_forces(s::SymbolicAWEModel) = s.get_spring_force(s.integrator)
function pos(s::SymbolicAWEModel)
    pos = s.get_pos(s.integrator)
    return [pos[:,i] for i in eachindex(pos[1,:])]
end    

function min_chord_length(s::SymbolicAWEModel)
    min_len = Inf
    for wing in s.vsm_wings
        le_pos = [wing.le_interp[i](wing.gamma_tip) for i in 1:3]
        te_pos = [wing.te_interp[i](wing.gamma_tip) for i in 1:3]
        min_len = min(norm(le_pos - te_pos), min_len)
    end
    return min_len
end

"""
    set_depower_steering!(s::SymbolicAWEModel, depower, steering) -> Nothing

Set kite depower and steering by adjusting tether lengths. Depower controls angle of attack,
steering controls left/right differential. Values are scaled by minimum chord length.
"""
function set_depower_steering!(s::SymbolicAWEModel, depower, steering)
    len = s.set_tether_length
    len .= tether_length(s)
    depower *= min_chord_length(s)
    steering *= min_chord_length(s)
    len[2] = 0.5 * (2*depower + 2*len[1] + steering)
    len[3] = 0.5 * (2*depower + 2*len[1] - steering)
    return nothing
end

"""
    set_v_wind_ground!(s::SymbolicAWEModel, v_wind_gnd=s.set.v_wind, upwind_dir=-π/2) -> Nothing

Set ground wind speed (m/s) and upwind direction (radians). Direction: 0=north, π/2=east, 
π=zouth, -π/2=west (default).
"""
function set_v_wind_ground!(s::SymbolicAWEModel, v_wind_gnd=s.set.v_wind, upwind_dir=-pi/2)
    s.set.v_wind = v_wind_gnd
    s.set.upwind_dir = rad2deg(upwind_dir)
    s.set_set(s.integrator, s.set)
    return nothing
end

function get_set_hash(set::Settings; 
        fields=[:segments, :model, :foil_file, :physical_model, :quasi_static, :winch_model]
    )
    h = zeros(UInt8, 1)
    for field in fields
        value = getfield(set, field)
        h = sha1(string((value, h)))
    end
    return h
end

function get_sys_struct_hash(sys_struct::SystemStructure)
    @unpack points, groups, segments, pulleys, tethers, winches, wings, transforms = sys_struct
    data_parts = []
    for point in points
        push!(data_parts, ("point", point.idx, point.wing_idx, Int(point.type)))
    end
    for segment in segments
        push!(data_parts, ("segment", segment.idx, segment.point_idxs, Int(segment.type)))
    end
    for group in groups
        push!(data_parts, ("group", group.idx, group.point_idxs, Int(group.type)))
    end
    for pulley in pulleys
        push!(data_parts, ("pulley", pulley.idx, pulley.segment_idxs, Int(pulley.type)))
    end
    for tether in tethers
        push!(data_parts, ("tether", tether.idx, tether.segment_idxs))
    end
    for winch in winches
        model_type = winch.model isa TorqueControlledMachine
        push!(data_parts, ("winch", winch.idx, model_type, winch.tether_idxs))
    end
    for wing in wings
        push!(data_parts, ("wing", wing.idx, wing.group_idxs))
    end
    for transform in transforms
        push!(data_parts, ("transform", transform.idx, transform.wing_idx, transform.rot_point_idx, 
                transform.base_point_idx, transform.base_transform_idx))
    end
    content = string(data_parts)
    return sha1(content)
end

# SPDX-License-Identifier: MIT

#= MIT License

Copyright (c) 2020, 2021, 2022, 2024 Uwe Fechner and Bart van de Lint

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

#= Models of a kite-power system in implicit form: residual = f(y, yd)

This model implements a 3D mass-spring system with reel-out. It uses six tether segments (the number can be
configured in the file data/settings.yaml). Two kite models are provided, the one point and the four point
kite model. The spring constant and the damping decrease with the segment length. The aerodynamic kite forces are
calculated, depending on reel-out speed, depower and steering settings. 

Scientific background: http://arxiv.org/abs/1406.6218 =#

module KiteModels

using PrecompileTools: @setup_workload, @compile_workload 
using Dierckx, Interpolations, Serialization, StaticArrays, LinearAlgebra, Parameters, NLsolve, Roots, ForwardDiff, DiffEqCallbacks,
      DocStringExtensions, OrdinaryDiffEqCore, OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, NonlinearSolve, FiniteDiff, DifferentiationInterface
using VortexStepMethod, ControlPlots
import Sundials
using Reexport, Pkg
@reexport using KitePodModels
@reexport using WinchModels
@reexport using AtmosphericModels
@reexport using Rotations
import Base.zero
import KiteUtils.calc_elevation
import KiteUtils.calc_heading
import KiteUtils.calc_course
import KiteUtils.SysState
import OrdinaryDiffEqCore.init
import OrdinaryDiffEqCore.step!
using ModelingToolkit, SymbolicIndexingInterface
using ModelingToolkit: t_nounits as t, D_nounits as D
using ADTypes: AutoFiniteDiff
using UnPack
import ModelingToolkit.SciMLBase: successful_retcode

export KPS3, KPS4, RamAirKite, KVec3, SimFloat, ProfileLaw, EXP, LOG, EXPLOG                     # constants and types
export calc_set_cl_cd!, copy_examples, copy_bin, update_sys_state!                            # helper functions
export clear!, find_steady_state!, residual!                                                  # low level workers
export init_sim!, init!, reinit!, next_step!, init_pos_vel, init_pos, model!                                 # high level workers
export pos_kite, calc_height, calc_elevation, calc_azimuth, calc_heading, calc_course, calc_orient_quat  # getters
export calc_azimuth_north, calc_azimuth_east
export winch_force, lift_drag, cl_cd, lift_over_drag, unstretched_length, tether_length, v_wind_kite     # getters
export calculate_rotational_inertia!
export kite_ref_frame, orient_euler, spring_forces, upwind_dir, copy_model_settings, menu2
import LinearAlgebra: norm

set_zero_subnormals(true)       # required to avoid drastic slow down on Intel CPUs when numbers become very small

# Constants
const G_EARTH = 9.81            # gravitational acceleration
const BRIDLE_DRAG = 1.1         # should probably be removed
const SYS_3L = "system_3l.yaml" # default system project for the 3L model

# Type definitions
"""
    const SimFloat = Float64

This type is used for all real variables, used in the Simulation. Possible alternatives: Float32, Double64, Dual
Other types than Float64 or Float32 do require support of Julia types by the solver. 
"""
const SimFloat = Float64

"""
   const KVec3    = MVector{3, SimFloat}

Basic 3-dimensional vector, stack allocated, mutable.
"""
const KVec3    = MVector{3, SimFloat}

"""
   const SVec3    = SVector{3, SimFloat}

Basic 3-dimensional vector, stack allocated, immutable.
"""
const SVec3    = SVector{3, SimFloat}  

"""
    abstract type AbstractKiteModel

All kite models must inherit from this type. All methods that are defined on this type must work
with all kite models. All exported methods must work on this type. 
"""
abstract type AbstractKiteModel end

"""
    const AKM = AbstractKiteModel

Short alias for the AbstractKiteModel. 
"""
const AKM = AbstractKiteModel

function __init__()
    if isdir(joinpath(pwd(), "data")) && isfile(joinpath(pwd(), "data", "system.yaml"))
        set_data_path(joinpath(pwd(), "data"))
    end
end

include("KPS4.jl") # include code, specific for the four point kite model
include("RamAirKite.jl") # include code, specific for the four point 3 line kite model
include("mtk_model.jl")
include("KPS3.jl") # include code, specific for the one point kite model
include("init.jl") # functions to calculate the initial state vector, the initial masses and initial springs
include("mtk_init.jl")
include("plotting.jl")

function menu2()
    Main.include("examples/menu2.jl")
end

# Calculate the lift and drag coefficient as a function of the angle of attack alpha.
function set_cl_cd!(s::AKM, alpha)
    angle =  rad2deg(alpha)
    s.alpha_2 = angle
    if angle > 180.0
        angle -= 360.0
    end
    if angle < -180.0
        angle += 360.0
    end
    s.param_cl = s.calc_cl(angle)
    s.param_cd = s.calc_cd(angle)
    nothing
end

# Calculate the angle of attack alpha from the apparend wind velocity vector
# v_app and the z unit vector of the kite reference frame.
function calc_alpha(v_app, vec_z)
    π/2.0 - acos(-(v_app ⋅ vec_z) / norm(v_app))
end

"""
    set_depower_steering!(s::AKM, depower, steering)

Setter for the depower and steering model inputs. 

Parameters:
- depower:   Relative depower,  must be between 0 .. 1.0
- steering:  Relative steering, must be between -1.0 .. 1.0.  

This function sets the variables s.depower, s.steering and s.alpha_depower. 

It takes the depower offset c0 and the dependency of the steering sensitivity from
the depower settings into account. The raw steering value is stored in s.kcu_steering.
"""
function set_depower_steering!(s::AKM, depower, steering)
    s.depower  = depower
    s.kcu_steering = steering
    s.alpha_depower = calc_alpha_depower(s.kcu, depower)
    s.steering = (steering - s.set.c0) / (1.0 + s.set.k_ds * (s.alpha_depower / deg2rad(s.set.alpha_d_max)))
    nothing
end


"""
    unstretched_length(s::AKM)

Getter for the unstretched tether reel-out length (at zero force).
"""
function unstretched_length(s::AKM) s.l_tether end


"""
    lift_drag(s::AKM)

Return a tuple of the scalar lift and drag forces. 

**Example:**  

    lift, drag = lift_drag(s)
"""
function lift_drag(s::AKM) return (norm(s.lift_force), norm(s.drag_force)) end

"""
    lift_over_drag(s::AKM)

Return the lift-over-drag ratio.
"""
function lift_over_drag(s::AKM)
    lift, drag = lift_drag(s)
    return lift / drag
end

"""
    v_wind_kite(s::AKM)

Return the vector of the wind speed at the height of the kite.
"""
function v_wind_kite(s::AKM) s.v_wind end

"""
    set_v_wind_ground!(s::AKM, height, v_wind_gnd=s.set.v_wind; upwind_dir=0.0)

Set the vector of the wind-velocity at the height of the kite. As parameter the height,
the ground wind speed [m/s] and the upwind direction [radians] are needed.
Is called by the function next_step!.
"""
function set_v_wind_ground!(s::AKM, height, v_wind_gnd=s.set.v_wind; upwind_dir=-pi/2)
    if height < 6.0
        height = 6.0
    end
    wind_dir = -upwind_dir - pi/2
    s.v_wind .= v_wind_gnd * calc_wind_factor(s.am, height) .* [cos(wind_dir), sin(wind_dir), 0]
    s.v_wind_gnd .= [v_wind_gnd * cos(wind_dir), v_wind_gnd * sin(wind_dir), 0.0]
    s.v_wind_tether .= s.v_wind_gnd * calc_wind_factor(s.am, height / 2.0)
    s.rho = calc_rho(s.am, height)
    nothing
end

function upwind_dir(s::AKM)
    upwind_dir(s.v_wind_gnd)
end
function upwind_dir(v_wind_gnd)
    if v_wind_gnd[1] == 0.0 && v_wind_gnd[2] == 0.0
        return NaN
    end
    wind_dir = atan(v_wind_gnd[2], v_wind_gnd[1])
    -(wind_dir + π/2)
end
@register_symbolic upwind_dir(v_wind_gnd)

"""
    tether_length(s::AKM)

Calculate and return the real, stretched tether length.
"""
function tether_length(s::AKM)
    length = 0.0
    for i in 1:s.set.segments
        length += norm(s.pos[i+1] - s.pos[i])
    end
    return length
end

"""
    orient_euler(s::AKM)

Calculate and return the orientation of the kite in euler angles (roll, pitch, yaw)
as SVector.
"""
function orient_euler(s::AKM; one_point=false)
    q = QuatRotation(calc_orient_quat(s; one_point))
    roll, pitch, yaw = quat2euler(q)
    SVector(roll, pitch, yaw)
end

function calc_orient_quat(s::AKM; viewer=false, one_point=false)
    if viewer
        x, _, z = kite_ref_frame(s)
        pos_kite_ = pos_kite(s)
        pos_before = pos_kite_ .+ z
    
        rotation = rot(pos_kite_, pos_before, -x)
    else
        x, y, z = kite_ref_frame(s; one_point) # in ENU reference
        x = enu2ned(x)
        y = enu2ned(y)
        z = enu2ned(z)
            
        # reference frame for the orientation: NED (north, east, down)
        ax = @SVector [1, 0, 0]
        ay = @SVector [0, 1, 0]
        az = @SVector [0, 0, 1]
        rotation = rot3d(ax, ay, az, x, y, z)
    end
    q = QuatRotation(rotation)
    return Rotations.params(q)
end

function calc_orient_quat_old(s::AKM)
    x, y, z = kite_ref_frame(s) # in ENU reference
        
    ax = [0, 1, 0] # in ENU reference frame this is pointing to the south
    ay = [1, 0, 0] # in ENU reference frame this is pointing to the west
    az = [0, 0, -1] # in ENU reference frame this is pointing down
    rotation = rot3d(ax, ay, az, x, y, z)
    q = QuatRotation(rotation)
    return Rotations.params(q)
end

"""
    calc_elevation(s::AKM)

Determine the elevation angle of the kite in radian.
"""
function calc_elevation(s::AKM)
    KiteUtils.calc_elevation(pos_kite(s))
end

"""
    calc_azimuth(s::AKM)

Determine the azimuth angle of the kite in wind reference frame in radian.
Positive anti-clockwise when seen from above.
"""
function calc_azimuth(s::AKM)
    azn = KiteUtils.azimuth_north(pos_kite(s))
    azn2azw(azn; upwind_dir = upwind_dir(s))
end

"""
    calc_azimuth_east(s::AKM)

Determine the azimuth_east angle of the kite in radian.
Positive clockwise when seen from above.
"""
function calc_azimuth_east(s::AKM)
    KiteUtils.azimuth_east(pos_kite(s))
end

"""
    calc_azimuth_north(s::AKM)

Determine the azimuth_north angle of the kite in radian.
Positive anti-clockwise when seen from above.
"""
function calc_azimuth_north(s::AKM)
    KiteUtils.azimuth_north(pos_kite(s))
end

"""
    calc_heading(s::AKM; upwind_dir_=upwind_dir(s), neg_azimuth=false, one_point=false)

Determine the heading angle of the kite in radian.
"""
function calc_heading(s::AKM; upwind_dir_=upwind_dir(s), neg_azimuth=false, one_point=false)
    orientation = orient_euler(s; one_point)
    elevation = calc_elevation(s)
    # use azimuth in wind reference frame
    if neg_azimuth 
        azimuth = -calc_azimuth(s)
    else
        azimuth = calc_azimuth(s)
    end
    calc_heading(orientation, elevation, azimuth; upwind_dir=upwind_dir_)
end

"""
    calc_course(s::AKM)

Determine the course angle of the kite in radian.
Undefined if the velocity of the kite is near zero.
"""
function calc_course(s::AKM, neg_azimuth=false)
    elevation = calc_elevation(s)
    if neg_azimuth 
        azimuth = -calc_azimuth(s)
    else    
        azimuth = calc_azimuth(s)
    end
    KiteUtils.calc_course(s.vel_kite, elevation, azimuth)
end

"""
    calculate_rotational_inertia!(s::AKM, include_kcu::Bool=true, around_kcu::Bool=false)

Calculate the rotational inertia (Ixx, Ixy, Ixz, Iyy, Iyz, Izz) for a kite model from settings. Modifies the kitemodel by initialising the masses.

Parameters:
- X: x-coordinates of the point masses.
- Y: y-coordinates of the point masses.
- Z: z-coordinates of the point masses.
- M: masses of the point masses.
- `include_kcu`: Include the kcu in the rotational intertia calculation?
- `around_kcu`: Uses the kcu as the rotation point.

Returns:  
The tuple  Ixx, Ixy, Ixz, Iyy, Iyz, Izz where:
- Ixx: rotational inertia around the x-axis.
- Ixy: rotational inertia around the xy-plane.
- Ixz: rotational inertia around the xz-plane.
- Iyy: rotational inertia around the y-axis.
- Iyz: rotational inertia around the yz-plane.
- Izz: rotational inertia around the z-axis.

"""
function calculate_rotational_inertia!(s::AKM, include_kcu::Bool=true, around_kcu::Bool=false)
    points = KiteUtils.get_particles(s.set.height_k, s.set.h_bridle, s.set.width, s.set.m_k, [0, 0, 0], [0, 0, -1], [10, 0, 0])
    
    pos_matrix = [points[begin+1] points[begin+2] points[begin+3] points[begin+4] points[begin+5]]
    X = pos_matrix[begin, :]
    Y = pos_matrix[begin+1, :]
    Z = pos_matrix[begin+2, :]

    masses = init_masses!(s)
    M = masses[s.set.segments+1:end]

    if !include_kcu
        X = X[begin+1:end]
        Y = Y[begin+1:end]
        Z = Z[begin+1:end]
        M = M[begin+1:end]
    end

    if around_kcu
        Ixx, Ixy, Ixz, Iyy, Iyz, Izz = KiteUtils.calculate_rotational_inertia(X, Y, Z, M, false, points[begin+1])
    else
        Ixx, Ixy, Ixz, Iyy, Iyz, Izz = KiteUtils.calculate_rotational_inertia(X, Y, Z, M)
    end

    Ixx, Ixy, Ixz, Iyy, Iyz, Izz
end

# mutable struct SysState{P}
#     "time since start of simulation in seconds"
#     time::Float64
#     "time needed for one simulation timestep"
#     t_sim::Float64
#     "state of system state control"
#     sys_state::Int16
#     "mechanical energy [Wh]"
#     e_mech::Float64
#     "orientation of the kite (quaternion, order w,x,y,z)"
#     orient::MVector{4, Float32}
#     "elevation angle in radians"
#     elevation::MyFloat
#     "azimuth angle in radians"
#     azimuth::MyFloat
#     "tether length [m]"
#     l_tether::MyFloat
#     "reel out velocity [m/s]"
#     v_reelout::MyFloat
#     "tether force [N]"
#     force::MyFloat
#     "depower settings [0..1]"
#     depower::MyFloat
#     "steering settings [-1..1]"
#     steering::MyFloat
#     "heading angle in radian"
#     heading::MyFloat
#     "course angle in radian"
#     course::MyFloat
#     "norm of apparent wind speed [m/s]"
#     v_app::MyFloat
#     "velocity vector of the kite"
#     vel_kite::MVector{3, MyFloat}
#     "vector of particle positions in x"
#     X::MVector{P, MyFloat}
#     "vector of particle positions in y"
#     Y::MVector{P, MyFloat}
#     "vector of particle positions in z"
#     Z::MVector{P, MyFloat}
#     var_01::MyFloat
#     var_02::MyFloat
#     ...
#     var_16::MyFloat
# end 

function update_sys_state!(ss::SysState, s::AKM, zoom=1.0)
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
    ss.force = winch_force(s)
    ss.heading = calc_heading(s)
    ss.course = calc_course(s)
    ss.v_app = norm(s.v_apparent)
    ss.l_tether = s.l_tether
    ss.v_reelout = s.v_reel_out
    ss.depower = s.depower
    ss.steering = s.steering/s.set.cs_4p
    ss.kcu_steering = s.kcu_steering/s.set.cs_4p
    ss.vel_kite .= s.vel_kite
    ss.t_sim = 0.0
    ss.AoA = deg2rad(s.alpha_2)
    if isa(s, KPS4)
        ss.alpha3 = deg2rad(s.alpha_3)
        ss.alpha4 = deg2rad(s.alpha_4)
        if isnothing(s.set_force)
            ss.set_force = NaN
        else
            ss.set_force = s.set_force
        end
        if isnothing(s.bearing)
            ss.bearing = NaN
        else
            ss.bearing = s.bearing
        end
        if isnothing(s.attractor)
            ss.attractor = [NaN, NaN]
        else
            ss.attractor = s.attractor
        end
    end
    ss.set_steering = s.kcu.set_steering
    if isnothing(s.set_torque)
        ss.set_torque = NaN
    else
        ss.set_torque = s.set_torque
    end
    if isnothing(s.sync_speed)
        ss.set_speed = NaN
    else
        ss.set_speed = s.sync_speed
    end
    ss.roll, ss.pitch, ss.yaw = orient_euler(s)
    cl, cd = cl_cd(s)
    ss.CL2 = cl
    ss.CD2 = cd
    ss.v_wind_gnd  .= s.v_wind_gnd
    ss.v_wind_200m .= s.v_wind_gnd * calc_wind_factor(s.am, 200.0)
    ss.v_wind_kite .= s.v_wind
    nothing
end

"""
    SysState(s::AKM, zoom=1.0)

Constructor for creating a SysState object from a kite model (KPS3 or KPS4).
The SysState object can be used either for logging or for displaying the
system state in a viewer. Optionally the position arrays can be zoomed
according to the requirements of the viewer.
"""
function SysState(s::AKM, zoom=1.0)
    ss = SysState{length(s.pos)}()
    update_sys_state!(ss, s, zoom)
    ss
end

function calc_pre_tension(s::AKM)
    forces = spring_forces(s)
    av_force = 0.0
    for i in 1:s.set.segments
        av_force += forces[i]
    end
    av_force /= s.set.segments
    res = av_force/s.set.c_spring
    if res < 0.0 res = 0.0 end
    if isnan(res) res = 0.0 end
    return res + 1.0
end

"""
    init_sim!(s::AKM; t_end=1.0, stiffness_factor=0.5, delta=0.001, upwind_dir=-pi/2, prn=false)

Initializes the integrator of the model.

Parameters:
- s:     an instance of an abstract kite model
- t_end: end time of the simulation; normally not needed
- stiffness_factor: factor applied to the tether stiffness during initialization
- delta: initial stretch of the tether during the steady state calculation
- upwind_dir: upwind direction in radians, the direction the wind is coming from. Zero is at north; 
              clockwise positive. Default: -pi/2, wind from west.
- prn: if set to true, print the detailed solver results

Returns:
An instance of a DAE integrator.
"""
function init_sim!(s::AKM; t_end=1.0, stiffness_factor=0.5, delta=0.0001, upwind_dir=-pi/2, prn=false)
    clear!(s)
    s.stiffness_factor = stiffness_factor
    
    try
        y0, yd0 = KiteModels.find_steady_state!(s; stiffness_factor, delta, upwind_dir, prn)
    catch e
        if e isa AssertionError
            println("ERROR: Failure to find initial steady state in find_steady_state! function!\n"*
                    "Try to increase the delta parameter or to decrease the initial_stiffness of the init_sim! function.")
            return nothing
        else
            rethrow(e)
        end
    end
    y0  = Vector{SimFloat}(y0)
    yd0 = Vector{SimFloat}(yd0)
    
    if s.set.solver=="IDA"
        solver  = Sundials.IDA(linear_solver=Symbol(s.set.linear_solver), max_order = s.set.max_order)
    elseif s.set.solver=="DImplicitEuler"
        solver  = DImplicitEuler(autodiff=AutoFiniteDiff())
    elseif s.set.solver=="DFBDF"
        solver  = DFBDF(autodiff=AutoFiniteDiff(), max_order=Val{s.set.max_order}())        
    else
        println("Error! Invalid solver in settings.yaml: $(s.set.solver)")
        return nothing
    end

    dt = 1/s.set.sample_freq
    tspan   = (0.0, dt) 
    abstol  = s.set.abs_tol # max error in m/s and m

    differential_vars = ones(Bool, length(y0))
    prob    = DAEProblem{true}(residual!, yd0, y0, tspan, s; differential_vars)
    integrator = OrdinaryDiffEqCore.init(prob, solver; abstol=abstol, reltol=s.set.rel_tol, save_everystep=false,
                                         initializealg=OrdinaryDiffEqCore.NoInit())
    if isa(s, KPS4)
        roll, pitch, yaw = orient_euler(s)
        s.pitch_rate = 0
        s.pitch = pitch
        set_initial_velocity!(s)
    end
    s.v_reel_out = s.set.v_reel_out
    return integrator
end


"""
    next_step!(s::AKM, integrator; set_speed = nothing, set_torque=nothing, set_force=nothing, bearing = nothing
               attractor=nothing, v_wind_gnd=s.set.v_wind, upwind_dir=-pi/2, dt=1/s.set.sample_freq)

Calculates the next simulation step. Either `set_speed` or `set_torque` must be provided.

Parameters:
- s:            an instance of an abstract kite model
- integrator:   an integrator instance as returned by the function [`init_sim!`](@ref)
- set_speed:         set value of reel out speed in m/s or nothing
- set_torque:   set value of the torque in Nm or nothing
- set_force:    set value of the force in N or nothing (only for logging, not used otherwise)
- bearing:      set value of heading/ course in radian or nothing (only for logging, not used otherwise)
- attractor:    the attractor coordinates [azimuth, elevation] in radian or nothing (only for logging)
- `v_wind_gnd`: wind speed at reference height in m/s
- `upwind_dir`: upwind direction in radians, the direction the wind is coming from. Zero is at north; 
                clockwise positive. Default: -pi/2, wind from west.
- dt:           time step in seconds

Returns:
The end time of the time step in seconds.
"""
function next_step!(s::AKM, integrator; set_speed = nothing, set_torque=nothing, set_force=nothing, bearing = nothing,
                    attractor=nothing, v_wind_gnd=s.set.v_wind, upwind_dir=-pi/2, dt=1/s.set.sample_freq)
    KitePodModels.on_timer(s.kcu)
    KiteModels.set_depower_steering!(s, get_depower(s.kcu), get_steering(s.kcu))
    s.sync_speed = set_speed
    s.set_torque = set_torque
    if isa(s, KPS4)
        s.set_force = set_force
        s.bearing = bearing
        s.attractor = attractor
    end
    s.t_0 = integrator.t
    set_v_wind_ground!(s, calc_height(s), v_wind_gnd; upwind_dir)
    s.iter = 0
    if s.set.solver == "IDA"
        Sundials.step!(integrator, dt, true)
    else
        OrdinaryDiffEqCore.step!(integrator, dt, true)
    end
    if s.stiffness_factor < 1.0
        s.stiffness_factor+=0.01
        if s.stiffness_factor > 1.0
            s.stiffness_factor = 1.0
        end
    end
    if isa(s, KPS4)
        roll, pitch, yaw = orient_euler(s)
        s.pitch_rate = (pitch - s.pitch) / dt
        s.pitch = pitch
    end
    integrator.t
end

"""
    copy_examples()

Copy all example scripts to the folder "examples"
(it will be created if it doesn't exist).
"""
function copy_examples()
    PATH = "examples"
    if ! isdir(PATH) 
        mkdir(PATH)
    end
    src_path = joinpath(dirname(pathof(@__MODULE__)), "..", PATH)
    copy_files("examples", readdir(src_path))
end

function copy_model_settings()
    files = ["settings.yaml", "MH82.dat", "polars.bin", "system.yaml", "settings_3l.yaml", 
             "system_3l.yaml"]
    dst_path = abspath(joinpath(pwd(), "data"))
    copy_files("data", files)
    set_data_path(joinpath(pwd(), "data"))
    println("Copied $(length(files)) files to $(dst_path) !")
end

function install_examples(add_packages=true)
    copy_examples()
    copy_settings()
    copy_model_settings()
    if add_packages
        Pkg.add("KiteUtils")
        Pkg.add("KitePodModels")
        Pkg.add("WinchModels")
        Pkg.add("ControlPlots")
        Pkg.add("LaTeXStrings")
        Pkg.add("StatsBase")
    end
end

function copy_files(relpath, files)
    if ! isdir(relpath) 
        mkdir(relpath)
    end
    src_path = joinpath(dirname(pathof(@__MODULE__)), "..", relpath)
    for file in files
        cp(joinpath(src_path, file), joinpath(relpath, file), force=true)
        chmod(joinpath(relpath, file), 0o774)
    end
    files
end

function create_bridle(se)
    # create the bridle
    bridle = KiteUtils.get_particles(se.height_k, se.h_bridle, se.width, se.m_k)
    return bridle
end
function bridle_length(se)
    # calculate the bridle length
    bridle = create_bridle(se)[2:end]
    len = norm(bridle[1] - bridle[2])
    len += norm(bridle[1] - bridle[4])
    len += norm(bridle[1] - bridle[5])
    len += norm(bridle[3] - bridle[2])
    len += norm(bridle[3] - bridle[4])
    len += norm(bridle[3] - bridle[5])
end


"""
    copy_bin()

Copy the scripts create_sys_image and run_julia to the folder "bin"
(it will be created if it doesn't exist).
"""
function copy_bin()
    PATH = "bin"
    if ! isdir(PATH) 
        mkdir(PATH)
    end
    src_path = joinpath(dirname(pathof(@__MODULE__)), "..", PATH)
    cp(joinpath(src_path, "create_sys_image"), joinpath(PATH, "create_sys_image"), force=true)
    cp(joinpath(src_path, "run_julia"), joinpath(PATH, "run_julia"), force=true)
    chmod(joinpath(PATH, "create_sys_image"), 0o774)
    chmod(joinpath(PATH, "run_julia"), 0o774)
    PATH = "test"
    if ! isdir(PATH) 
        mkdir(PATH)
    end
    src_path = joinpath(dirname(pathof(@__MODULE__)), "..", PATH)
    cp(joinpath(src_path, "create_sys_image2.jl"), joinpath(PATH, "create_sys_image.jl"), force=true)
    cp(joinpath(src_path, "test_for_precompile.jl"), joinpath(PATH, "test_for_precompile.jl"), force=true)
    cp(joinpath(src_path, "update_packages.jl"), joinpath(PATH, "update_packages.jl"), force=true)
    chmod(joinpath(PATH, "create_sys_image.jl"), 0o664)
    chmod(joinpath(PATH, "test_for_precompile.jl"), 0o664)
    chmod(joinpath(PATH, "update_packages.jl"), 0o664)
end

# @setup_workload begin
#     # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
#     # precompile file and potentially make loading faster.
#     # list = [OtherType("hello"), OtherType("world!")]
#     path = dirname(pathof(@__MODULE__))
#     set_data_path(joinpath(path, "..", "data"))

#     set = se("system.yaml")
#     set.kcu_diameter = 0
#     kps4_::KPS4 = KPS4(KCU(set))
#     kps3_::KPS3 = KPS3(KCU(se("system.yaml")))
#     if ! haskey(ENV, "NO_MTK")    
#         kpsq_::RamAirKite = RamAirKite(KCU(se(SYS_3L)))
#     end
#     @assert ! isnothing(kps4_.wm)
#     @compile_workload begin
#         # all calls in this block will be precompiled, regardless of whether
#         # they belong to your package or not (on Julia 1.8 and higher)
#         integrator = KiteModels.init_sim!(kps3_; stiffness_factor=0.035, prn=false)
#         integrator = KiteModels.init_sim!(kps4_; delta=0.03, stiffness_factor=0.05, prn=false) 
#         if ! haskey(ENV, "NO_MTK")
#             integrator = KiteModels.init_sim!(kpsq_)
#         end
#         nothing
#     end
# end
end
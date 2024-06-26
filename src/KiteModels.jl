# SPDX-License-Identifier: MIT

#= MIT License

Copyright (c) 2020, 2021, 2022 Uwe Fechner

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
using Dierckx, StaticArrays, Rotations, LinearAlgebra, Parameters, NLsolve, DocStringExtensions, OrdinaryDiffEq
import Sundials
using Reexport
@reexport using KitePodModels
@reexport using WinchModels
@reexport using AtmosphericModels
import Base.zero
import KiteUtils.calc_elevation
import KiteUtils.calc_azimuth
import KiteUtils.calc_heading
import KiteUtils.calc_course
import KiteUtils.SysState
# import Sundials.init
# import Sundials.step!
import OrdinaryDiffEq.init
import OrdinaryDiffEq.step!

export KPS3, KPS4, KPS4_3L, KVec3, SimFloat, ProfileLaw, EXP, LOG, EXPLOG                              # constants and types
export calc_set_cl_cd!, copy_examples, copy_bin, update_sys_state!                            # helper functions
export clear!, find_steady_state!, residual!                                                  # low level workers
export init_sim!, next_step!                                                                  # high level workers
export pos_kite, calc_height, calc_elevation, calc_azimuth, calc_heading, calc_course         # getters
export winch_force, lift_drag, lift_over_drag, unstretched_length, tether_length, v_wind_kite # getters
export kite_ref_frame, orient_euler, spring_forces
import LinearAlgebra: norm

set_zero_subnormals(true)       # required to avoid drastic slow down on Intel CPUs when numbers become very small

# Constants
const G_EARTH = 9.81            # gravitational acceleration
const BRIDLE_DRAG = 1.1         # should probably be removed

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

# the following two definitions speed up the function residual! from 940ns to 540ns
# disadvantage: changing the cl and cd curves requires a restart of the program     
const calc_cl = Spline1D(se().alpha_cl, se().cl_list)
const calc_cd = Spline1D(se().alpha_cd, se().cd_list)  

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
include("KPS4_3L.jl") # include code, specific for the four point kite model
include("KPS3.jl") # include code, specific for the one point kite model
include("init.jl") # functions to calculate the inital state vector, the inital masses and initial springs

# Calculate the lift and drag coefficient as a function of the angle of attack alpha.
function set_cl_cd!(s::AKM, alpha)   
    angle =  alpha * 180.0 / π
    if angle > 180.0
        angle -= 360.0
    end
    if angle < -180.0
        angle += 360.0
    end
    s.param_cl = calc_cl(angle)
    s.param_cd = calc_cd(angle)
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
the depower settings into account.
"""
function set_depower_steering!(s::AKM, depower, steering)
    s.depower  = depower
    s.alpha_depower = calc_alpha_depower(s.kcu, depower)
    s.steering = (steering - s.set.c0) / (1.0 + s.set.k_ds * (s.alpha_depower / deg2rad(s.set.alpha_d_max)))
    nothing
end


"""
    set_v_reel_out!(s::AKM, v_reel_out, t_0, period_time = 1.0 / s.set.sample_freq)

Setter for the reel-out speed. Must be called on every timestep (before each simulation).
It also updates the tether length, therefore it must be called even if `v_reel_out` has
not changed.

- t_0 the start time of the next timestep relative to the start of the simulation [s]
"""
function set_v_reel_out!(s::AKM, v_reel_out, t_0, period_time = 1.0 / s.set.sample_freq)
    s.sync_speed = v_reel_out
    s.last_v_reel_out = s.v_reel_out
    s.t_0 = t_0
end

"""
    unstretched_length(s::AKM)

Getter for the unstretched tether reel-out lenght (at zero force).
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
    set_v_wind_ground!(s::AKM, height, v_wind_gnd=s.set.v_wind, wind_dir=0.0)

Set the vector of the wind-velocity at the height of the kite. As parameter the height,
the ground wind speed [m/s] and the wind direction [radians] are needed.
Must be called every at each timestep.
"""
function set_v_wind_ground!(s::AKM, height, v_wind_gnd=s.set.v_wind, wind_dir=0.0)
    if height < 6.0
        height = 6.0
    end
    s.v_wind .= v_wind_gnd * calc_wind_factor(s.am, height) .* [cos(wind_dir), sin(wind_dir), 0]
    s.v_wind_gnd .= [v_wind_gnd * cos(wind_dir), v_wind_gnd * sin(wind_dir), 0.0]
    s.v_wind_tether .= v_wind_gnd * calc_wind_factor(s.am, height / 2.0) .* [cos(wind_dir), sin(wind_dir), 0]
    s.rho = calc_rho(s.am, height)
    nothing
end

"""
    tether_length(s::AKM)

Calculate and return the real, stretched tether lenght.
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
function orient_euler(s::AKM)
    x, y, z = kite_ref_frame(s)
    roll = atan(y[3], z[3]) - π/2
    if roll < -π/2
       roll += 2π
    end
    pitch = asin(-x[3])
    yaw = -atan(x[2], x[1]) - π/2
    if yaw < -π/2
        yaw += 2π
    end
    SVector(roll, pitch, yaw)
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

Determine the azimuth angle of the kite in radian.
"""
function calc_azimuth(s::AKM)
    KiteUtils.azimuth_east(pos_kite(s))
end

"""
    calc_heading(s::AKM)

Determine the heading angle of the kite in radian.
"""
function calc_heading(s::AKM)
    orientation = orient_euler(s)
    elevation = calc_elevation(s)
    azimuth = calc_azimuth(s)
    KiteUtils.calc_heading(orientation, elevation, azimuth)
end

"""
    calc_course(s::AKM)

Determine the course angle of the kite in radian.
Undefined if the velocity of the kite is near zero.
"""
function calc_course(s::AKM)
    elevation = calc_elevation(s)
    azimuth = calc_azimuth(s)
    KiteUtils.calc_course(s.vel_kite, elevation, azimuth)
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
#     var_03::MyFloat
#     var_04::MyFloat
#     var_05::MyFloat
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
    x, y, z = kite_ref_frame(s)
    pos_kite_ = pos_kite(s)
    pos_before = pos_kite_ + z
   
    rotation = rot(pos_kite_, pos_before, -x)
    q = QuatRotation(rotation)
    ss.orient .= Rotations.params(q)
    ss.elevation = calc_elevation(s)
    ss.azimuth = calc_azimuth(s)
    ss.force = winch_force(s)
    ss.heading = calc_heading(s)
    ss.course = calc_course(s)
    ss.v_app = norm(s.v_apparent)
    ss.l_tether = s.l_tether
    ss.v_reelout = s.v_reel_out
    ss.depower = s.depower
    ss.steering = s.steering
    ss.vel_kite .= s.vel_kite
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
    
    x, y, z = kite_ref_frame(s)
    pos_kite_ = pos_kite(s)
    pos_before = pos_kite_ + z
   
    rotation = rot(pos_kite_, pos_before, -x)
    q = QuatRotation(rotation)
    orient = MVector{4, Float32}(Rotations.params(q))

    elevation = calc_elevation(s)
    azimuth = calc_azimuth(s)
    force = winch_force(s)
    heading = calc_heading(s)
    course = calc_course(s)
    v_app_norm = norm(s.v_apparent)
    t_sim = 0
    KiteUtils.SysState{P}(s.t_0, t_sim, 0, 0, orient, elevation, azimuth, s.l_tether, s.v_reel_out, force, s.depower, s.steering, 
                          heading, course, v_app_norm, s.vel_kite, X, Y, Z, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
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
    init_sim!(s; t_end=1.0, stiffness_factor=0.035, prn=false)

Initialises the integrator of the model.

Parameters:
- s:     an instance of an abstract kite model
- t_end: end time of the simulation; normally not needed
- stiffness_factor: factor applied to the tether stiffness during initialisation
- prn: if set to true, print the detailed solver results

Returns:
An instance of a DAE integrator.
"""
function init_sim!(s::AKM; t_end=1.0, stiffness_factor=0.035, prn=false)
    clear!(s)
    s.stiffness_factor = stiffness_factor
    y0, yd0 = KiteModels.find_steady_state!(s; stiffness_factor=stiffness_factor, prn=prn)
    y0  = Vector{Float64}(y0)
    yd0 = Vector{Float64}(yd0)
    if s.set.solver=="IDA"
        solver  = Sundials.IDA(linear_solver=Symbol(s.set.linear_solver), max_order = s.set.max_order)
    elseif s.set.solver=="DImplicitEuler"
        solver  = DImplicitEuler(autodiff=false)
    elseif s.set.solver=="DFBDF"
        solver  = DFBDF(autodiff=false, max_order=Val{s.set.max_order}())
    else
        println("Error! Invalid solver in settings.yaml: $(s.set.solver)")
        return nothing
    end
    dt = 1/s.set.sample_freq
    tspan   = (0.0, dt) 
    abstol  = s.set.abs_tol # max error in m/s and m
    differential_vars = ones(Bool, length(y0))
    prob    = DAEProblem{true}(residual!, yd0, y0, tspan, s; differential_vars)
    integrator = OrdinaryDiffEq.init(prob, solver; abstol=abstol, reltol= s.set.rel_tol, save_everystep=false)
end

"""
    next_step!(s::AKM, integrator; v_ro = 0.0, v_wind_gnd=s.set.v_wind, wind_dir=0.0, dt=1/s.set.sample_freq)

Calculates the next simulation step.

Parameters:
- s:            an instance of an abstract kite model
- integrator:   an integrator instance as returned by the function [`init_sim!`](@ref)
- v_ro:         set value of reel out speed in m/s
- `v_wind_gnd`: wind speed at reference height in m/s
- wind_dir:     wind direction in radians
- dt:           time step in seconds

Only the first two parameters are required.

Returns:
The end time of the time step in seconds.
"""
function next_step!(s::AKM, integrator; v_ro = 0.0, v_wind_gnd=s.set.v_wind, wind_dir=0.0, dt=1/s.set.sample_freq)
    KitePodModels.on_timer(s.kcu)
    KiteModels.set_depower_steering!(s, get_depower(s.kcu), get_steering(s.kcu))
    set_v_reel_out!(s, v_ro, integrator.t)
    set_v_wind_ground!(s, calc_height(s), v_wind_gnd, wind_dir)
    if s.set.solver == "IDA"
        Sundials.step!(integrator, dt, true)
    else
        OrdinaryDiffEq.step!(integrator, dt, true)
    end
    if s.stiffness_factor < 1.0
        s.stiffness_factor+=0.01
        if s.stiffness_factor > 1.0
            s.stiffness_factor = 1.0
        end
    end
    integrator.t
end

"""
    copy_examples()

Copy the example scripts to the folder "examples"
(it will be created if it doesn't exist).
"""
function copy_examples()
    PATH = "examples"
    if ! isdir(PATH) 
        mkdir(PATH)
    end
    src_path = joinpath(dirname(pathof(@__MODULE__)), "..", PATH)
    cp(joinpath(src_path, "compare_kps3_kps4.jl"), joinpath(PATH, "compare_kps3_kps4.jl"), force=true)
    cp(joinpath(src_path, "simulate.jl"), joinpath(PATH, "simulate.jl"), force=true)
    cp(joinpath(src_path, "reel_out_1p.jl"), joinpath(PATH, "reel_out_1p.jl"), force=true)
    cp(joinpath(src_path, "reel_out_4p.jl"), joinpath(PATH, "reel_out_4p.jl"), force=true)
    chmod(joinpath(PATH, "compare_kps3_kps4.jl"), 0o664)
    chmod(joinpath(PATH, "simulate.jl"), 0o664)
    chmod(joinpath(PATH, "reel_out_1p.jl"), 0o664)
    chmod(joinpath(PATH, "reel_out_4p.jl"), 0o664)
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
    cp(joinpath(src_path, "create_sys_image2"), joinpath(PATH, "create_sys_image"), force=true)
    cp(joinpath(src_path, "run_julia2"), joinpath(PATH, "run_julia"), force=true)
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

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    # list = [OtherType("hello"), OtherType("world!")]
    set_data_path()
    kps4_::KPS4 = KPS4(KCU(se()))
    kps3_::KPS3 = KPS3(KCU(se()))
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        integrator = KiteModels.init_sim!(kps3_; stiffness_factor=0.035, prn=false)
        integrator = KiteModels.init_sim!(kps4_; stiffness_factor=0.035, prn=false)       
        nothing
    end
end
end
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

#= Model of a kite-power system in implicit form: residual = f(y, yd)

This model implements a 3D mass-spring system with reel-out. It uses five tether segments (the number can be
configured in the file data/settings.yaml). The kite is modelled as additional mass at the end of the tether.
The spring constant and the damping decrease with the segment length. The aerodynamic kite forces are
calculated, depending on reel-out speed, depower and steering settings. 

Scientific background: http://arxiv.org/abs/1406.6218 =#

module KiteModels

using Dierckx, StaticArrays, LinearAlgebra, Parameters, NLsolve, DocStringExtensions
using KiteUtils, KitePodModels
import Base.zero

export KPS3, KPS4, KVec3, SimFloat, ProfileLaw, EXP, LOG, EXPLOG                                  # constants and types
export calc_rho, calc_wind_factor, calc_drag, calc_set_cl_cd, clear, find_steady_state, residual! # environment and helper functions
export calc_height                                                                                # helper functions
export set_v_reel_out, set_depower_steering                                                       # setters
export winch_force, lift_drag, lift_over_drag, unstretched_length, tether_length, v_wind_kite     # getters
export spring_forces

set_zero_subnormals(true)           # required to avoid drastic slow down on Intel CPUs when numbers become very small

# Constants
const G_EARTH = 9.81                # gravitational acceleration
const BRIDLE_DRAG = 1.1             # should probably be removed

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

include("KPS4.jl") # include code, specific for the four point kite model
include("KPS3.jl") # include code, specific for the one point kite model

"""
    calc_rho(s, height)

Calculate the air densisity as function of height.
"""
@inline function calc_rho(s::AKM, height) s.set.rho_0 * exp(-height / 8550.0) end

"""
    ProfileLaw

Enumeration to describe the wind profile low that is used.
"""
@enum ProfileLaw EXP=1 LOG=2 EXPLOG=3

"""
    calc_wind_factor(s, height, profile_law=s.set.profile_law)

Calculate the relative wind speed at a given height and reference height.
"""
@inline function calc_wind_factor(s, height, profile_law=s.set.profile_law)
    if typeof(profile_law) != ProfileLaw
        profile_law = ProfileLaw(profile_law)
    end
    if height < s.set.h_ref
        height = s.set.h_ref
    end
    if profile_law == EXP
        return exp(s.set.alpha * log(height/s.set.h_ref))
    elseif profile_law == LOG
        return log(height / s.set.z0) / log(s.set.h_ref / s.set.z0)
    else
        K = 1.0
        log1 = log(height / s.set.z0) / log(s.set.h_ref / s.set.z0)
        exp1 = exp(s.set.alpha * log(height/s.set.h_ref))
        return log1 +  K * (log1 - exp1)
    end
end

# Calculate the lift and drag coefficient as a function of the angle of attack alpha.
function set_cl_cd(s::AKM, alpha)   
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
    π/2.0 - acos(-dot(v_app, vec_z) / norm(v_app))
end

"""
    calc_set_cl_cd(s, vec_c, v_app)

Calculate the lift over drag ratio as a function of the direction vector of the last tether
segment, the current depower setting and the apparent wind speed.
Set the calculated CL and CD values in the struct s. 
"""
function calc_set_cl_cd(s::AKM, vec_c, v_app)
    s.vec_z .= normalize(vec_c)
    alpha = calc_alpha(v_app, s.vec_z) - s.alpha_depower
    set_cl_cd(s, alpha)
end

"""
    set_depower_steering(s::KPS3, depower, steering)

Setter for the depower and steering model inputs. 
- valid range for steering: -1.0 .. 1.0.  
- valid range for depower: 0 .. 1.0

This function sets the variables s.depower, s.steering and s.alpha_depower. 

It takes the depower offset c0 and the dependency of the steering sensitivity from
the depower settings into account.
"""
function set_depower_steering(s::AKM, depower, steering)
    s.depower  = depower
    s.alpha_depower = calc_alpha_depower(s.kcu, depower)
    s.steering = (steering - s.set.c0) / (1.0 + s.set.k_ds * (s.alpha_depower / deg2rad(s.set.alpha_d_max)))
    nothing
end


"""
    set_v_reel_out(s::AKM, v_reel_out, t_0, period_time = 1.0 / s.set.sample_freq)

Setter for the reel-out speed. Must be called on every timestep (before each simulation).
It also updates the tether length, therefore it must be called even if `v_reel_out` has
not changed.

- t_0 the start time of the next timestep relative to the start of the simulation [s]
"""
function set_v_reel_out(s::AKM, v_reel_out, t_0, period_time = 1.0 / s.set.sample_freq)
    s.l_tether += 0.5 * (v_reel_out + s.last_v_reel_out) * period_time
    s.last_v_reel_out = s.v_reel_out
    s.v_reel_out = v_reel_out
    s.t_0 = t_0
end

"""
    unstretched_length(s::AKM)

Getter for the unstretched tether reel-out lenght (at zero force).
"""
function unstretched_length(s::AKM) s.l_tether end

"""
    winch_force(s::AKM)

Return the absolute value of the force at the winch as calculated during the last timestep. 
"""
function winch_force(s::AKM) norm(s.last_force) end

"""
    spring_forces(s::AKM)

Return an array of the scalar spring forces of all tether segements.
"""
function spring_forces(s::AKM)
    forces = zeros(SimFloat, s.set.segments)
    for i in 1:s.set.segments
        forces[i] =  s.c_spring * (norm(s.pos[i+1] - s.pos[i]) - s.length)
    end
    forces
end

"""
    lift_drag(s::AKM)

Return a tuple of the scalar lift and drag forces. 

**Example:**  

    lift, drag = lift_drag(kps)
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
    set_v_wind_ground(s::AKM, height, v_wind_gnd=s.set.v_wind, wind_dir=0.0)

Set the vector of the wind-velocity at the height of the kite. As parameter the height,
the ground wind speed and the wind direction are needed.
Must be called every at each timestep.
"""
function set_v_wind_ground(s::AKM, height, v_wind_gnd=s.set.v_wind, wind_dir=0.0)
    if height < 6.0
        height = 6.0
    end
    s.v_wind .= v_wind_gnd * calc_wind_factor(s, height) .* [cos(wind_dir), sin(wind_dir), 0]
    s.v_wind_gnd .= [v_wind_gnd * cos(wind_dir), v_wind_gnd * sin(wind_dir), 0.0]
    s.v_wind_tether .= v_wind_gnd * calc_wind_factor(s, height / 2.0) .* [cos(wind_dir), sin(wind_dir), 0]
    s.rho = calc_rho(s, height)
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

precompile(find_steady_state, (KPS3{SimFloat, KVec3, 7},))   

end
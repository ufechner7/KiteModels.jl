using DiscretePIDs
# low level winch controller; this code will be moved to WinchControllers.jl in the future

mutable struct WinchSpeedController
    kp::Float64
    ki::Float64
    pid::DiscretePID
end
function WinchSpeedController(;kp=20.0, ki=5.0, dt)
    pid = DiscretePID(;K=kp, Ti=kp/ki, Ts=dt)
    WinchSpeedController(kp, ki, pid)
end
"""
    calc_set_torque(set::Settings, wcs::WinchSpeedController, v_set, v_reelout, force)

Calculate the set torque for the winch as defined in settings.yaml file.

# Arguments
- `set::Settings`: the settings struct
- `wcs::WinchSpeedController`: the winch speed controller
- `v_set`: the setpoint of the reelout speed
- `v_reelout`: the current reelout speed
- `force`: the current tether force
"""
function calc_set_torque(set::Settings, wcs::WinchSpeedController, v_set, v_reelout, force)
    # calculate the set force
    set_force = DiscretePIDs.calculate_control!(wcs.pid, v_set, v_reelout)
    err = set_force - force
    # calculate the set torque
    r = set.drum_radius
    n = set.gear_ratio
    set_torque = -r/n * (0.0*set_force-err)
end
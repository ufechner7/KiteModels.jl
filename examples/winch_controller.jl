using DiscretePIDs
# low level winch controller; this code will be moved to WinchControllers.jl in the future

mutable struct WinchSpeedController
    kp::Float64
    ki::Float64
    pid::DiscretePID
end
function WinchSpeedController(;kp=10.0, ki=0.0, dt)
    pid = DiscretePID(;K=kp, Ti=kp/ki, dt)
    WinchSpeedController(kp, ki, pid)
end
"""
    calc_set_torque(set::Settings, v_set, v_reelout, force)

Calculate the set torque for the winch as defined in settings.yaml file.

# Arguments
- `set::Settings`: the settings struct
- `v_set`: the setpoint of the reelout speed
- `v_reelout`: the current reelout speed
- `force`: the current tether force
"""
function calc_set_torque(set::Settings, wcs::WinchSpeedController, v_set, v_reelout, force)
    # calculate the speed error
    err = v_set - v_reelout
    # calculate the set force
    set_force = DiscretePIDs.update(wcs.pid, err)
    err = set_force - force
    # calculate the set torque
    r = set.drum_radius
    n = set.gear_ratio
    set_torque = -r/n * err
end
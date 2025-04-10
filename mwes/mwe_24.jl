"""
State estimation

function f(u)
    wind_vel, wind_dir, distance_acc = u
    init(wind_vel, wind_dir, distance_acc, point_acc = 0, twist_acc = 0)
    return state .- measured # where the chosen state and measured are of size u, and these are coupled to u
end
"""
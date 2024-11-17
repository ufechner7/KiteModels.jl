using Revise, KiteModels, ModelingToolkit, NonlinearSolve

"""
No need for lift information to be able to guess tether amount of slack
We have the actual tether length, torque, acc and tether force.
    Optimize flap_angle so that flap_acc is zero.
    Optimize the distance of the kite from ground station and the amount of 
    tether curve for each tether so that tether acceleration is closest to the 
    measured tether acceleration.

elevation ~ measured_elevation
azimuth ~ measured_azimuth
tether_length ~ measured_tether_length
tether_force ~ measured_tether_force


Update the init_pos_vel_acc function so that it takes in measurements and optimizing variables and only creates pos. Then
    use the found pos in the existing model to calculate accelerations. Optimize for least acceleration.
"""


set = deepcopy(load_settings("system_3l.yaml"))
# set.elevation = 71
dt = 0.05
total_time = 10.0

steps = Int(round(total_time / dt))
logger = Logger(3*set.segments + 6, steps)

if !@isdefined(s) s::KPS4_3L = KPS4_3L(KCU(set)) end
s.set = update_settings()
s.set.abs_tol = 0.0006
s.set.rel_tol = 0.001
s.set.l_tether = 21.0
# s.set.damping = 473
s.set.elevation = 87
s.torque_control = true

ns = nonlin_model!(s, s.pos, s.vel)

guesses = Dict()
[guesses[ns.segment[i, j]] = 1 for i in 1:3 for j in 1:s.num_A]
[guesses[ns.γ[i, j]] = 1 for i in 1:3 for j in 1:s.set.segments-1]
[guesses[ns.mass_tether_particle[i]] = 1 for i in 1:3]

prob = NonlinearProblem(nsys, guesses, []; check_length=false)
sol = solve(prob, NewtonRaphson())

# guesses = defaults(ns)

# function get_guesses(s::KPS4_3L)
#     for i in eachindex(s.springs)
#         global guesses
#         p1 = s.springs[i].p1
#         p2 = s.springs[i].p2
#         for j in 1:3
#             guesses[ns.segment[j, i]] = pos[p1][j] - pos[p2][j]
#         end
#     end
#     for (j, i) in enumerate(s.set.segments*3-2:s.set.segments*3)
#         guesses[ns.l_0[i]] = s.tether_lengths[j] / s.set.segments
#     end
#     for i in 1:s.num_A
#         for j in 1:3
#             guesses[ns.v_apparent[j, i]] = s.v_wind[j]
#         end
#     end
#     P_c = 0.5 * (pos[s.num_D] + pos[s.num_C])
#     [guesses[ns.P_c[i]] = P_c[i] for i in 1:3]

#     n = s.set.aero_surfaces
#     α_0         = π/2 - s.set.width/2/s.set.radius
#     α_middle    = π/2
#     dα          = (α_middle - α_0) / n
#     KiteModels.calc_kite_ref_frame!(s, pos[s.num_E], pos[s.num_C], pos[s.num_D])
#     for i in 1:2n
#         if i <= n
#             α = α_0 + -dα/2 + i * dα
#         else
#             α = pi - (α_0 + -dα/2 + (i-n) * dα)
#         end
#         guesses[ns.seg_flap_angle[i]] = deg2rad(s.set.alpha_zero)
#         [guesses[ns.F[j, i]] = pos[s.num_E][j] + s.e_z[j] * (-s.set.bridle_center_distance + s.set.radius) + 
#             s.e_y[j] * cos(α) * s.set.radius - s.e_z[j] * sin(α) * s.set.radius for j in 1:3]
#     end
#     return guesses
# end
# guesses = get_guesses(s)
# prob = NonlinearProblem(ns, guesses, [ns.set_values[i] => set_values[i] for i in 1:3])
# @time sol = solve(prob, NewtonRaphson(); abstol=1e-3)

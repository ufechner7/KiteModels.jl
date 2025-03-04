# ==================== mtk model functions ================================================
# Implementation of the three-line model using ModellingToolkit.jl

function calc_acc_speed(motor::AsyncMachine, tether_vel, norm_, set_speed)
    calc_acceleration(motor, tether_vel, norm_; set_speed, set_torque=nothing, use_brake=false) # TODO: add brake setting
end
@register_symbolic calc_acc_speed(motor::AsyncMachine, tether_vel, norm_, set_speed)

function calc_acc_torque(motor::TorqueControlledMachine, tether_vel, norm_, set_torque)
    calc_acceleration(motor, tether_vel, norm_; set_speed=nothing, set_torque, use_brake=false)
end
@register_symbolic calc_acc_torque(motor::TorqueControlledMachine, tether_vel, norm_, set_torque)

function sym_interp(interp::Function, aoa, trailing_edge_angle)
    return interp(rad2deg(aoa), rad2deg(trailing_edge_angle-aoa)) # TODO: register callable struct https://docs.sciml.ai/Symbolics/dev/manual/functions/#Symbolics.@register_array_symbolic
end
@register_symbolic sym_interp(interp::Function, aoa, trailing_edge_angle)

function normalize(vec)
    return vec / norm(vec)
end
@register_symbolic normalize(vec)

function mean(vec)
    return sum(vec) / length(vec)
end

function convert_pos_vel(s::KPSQ, pos_, vel_)
    pos = Array{Union{Nothing, SimFloat}}(nothing, 3, s.i_A)
    vel = Array{Union{Nothing, SimFloat}}(nothing, 3, s.i_A)
    [pos[:,i] .= pos_[i] for i in 1:s.i_A-1]
    [vel[:,i] .= vel_[i] for i in 1:s.i_A-1]
    [pos[:,i] .= pos_[i] for i in s.i_B+1:s.i_A]
    [vel[:,i] .= vel_[i] for i in s.i_B+1:s.i_A]
    return pos, vel
end

function rotate_by_quaternion(q, v)
    p = [0, v...]
    return quaternion_multiply(q, quaternion_multiply(p, quaternion_conjugate(q)))[2:4]
end

function quaternion_conjugate(q)
    return [q[1], -q[2], -q[3], -q[4]]
end

function quaternion_multiply(q1, q2)
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1*w2 - x1*x2 - y1*y2 - z1*z2
    x = w1*x2 + x1*w2 + y1*z2 - z1*y2
    y = w1*y2 - x1*z2 + y1*w2 + z1*x2
    z = w1*z2 + x1*y2 - y1*x2 + z1*w2
    return [w, x, y, z]
end

function quaternion_to_rotation_matrix(q)
    w, x, y, z = q[1], q[2], q[3], q[4]
    
    return [
        1 - 2*(y*y + z*z)  2*(x*y - z*w)      2*(x*z + y*w);
        2*(x*y + z*w)      1 - 2*(x*x + z*z)  2*(y*z - x*w);
        2*(x*z - y*w)      2*(y*z + x*w)      1 - 2*(x*x + y*y)
    ]
end

function rotation_matrix_to_quaternion(R)
    tr_ = tr(R)
    
    if tr_ > 0
        S = sqrt(tr_ + 1.0) * 2
        w = 0.25 * S
        x = (R[3,2] - R[2,3]) / S
        y = (R[1,3] - R[3,1]) / S
        z = (R[2,1] - R[1,2]) / S
    elseif (R[1,1] > R[2,2]) && (R[1,1] > R[3,3])
        S = sqrt(1.0 + R[1,1] - R[2,2] - R[3,3]) * 2
        w = (R[3,2] - R[2,3]) / S
        x = 0.25 * S
        y = (R[1,2] + R[2,1]) / S
        z = (R[1,3] + R[3,1]) / S
    elseif R[2,2] > R[3,3]
        S = sqrt(1.0 + R[2,2] - R[1,1] - R[3,3]) * 2
        w = (R[1,3] - R[3,1]) / S
        x = (R[1,2] + R[2,1]) / S
        y = 0.25 * S
        z = (R[2,3] + R[3,2]) / S
    else
        S = sqrt(1.0 + R[3,3] - R[1,1] - R[2,2]) * 2
        w = (R[2,1] - R[1,2]) / S
        x = (R[1,3] + R[3,1]) / S
        y = (R[2,3] + R[3,2]) / S
        z = 0.25 * S
    end
    
    return [w, x, y, z]
end

# """ 
# Calculate the drag force and spring force of the tether segment, defined by the parameters pos1, pos2, vel1 and vel2
# and distribute it equally on the two particles, that are attached to the segment.
# The result is stored in the array s.forces. 
# """
# function calc_particle_forces!(s::KPSQ, eqs, force_eqs, force, p1, p2, pos1, pos2, vel1, vel2, length, c_spring, 
#     damping, rho, i, l_0, k, c, segment, rel_vel, av_vel, norm1, unit_vector, k2, c1, c2, spring_vel,
#             spring_force, v_apparent, v_wind_tether, area, v_app_perp, half_drag_force)
#     d_tether = s.set.d_tether/1000.0
#     eqs = [
#         eqs
#         l_0 ~ length[(i-1) % 3 + 1] # Unstressed length
#         k   ~ c_spring[(i-1) % 3 + 1] # Spring constant
#         c   ~ damping[(i-1) % 3 + 1] # Damping coefficient    
#         segment     .~ pos1 - pos2
#         rel_vel     .~ vel1 - vel2
#         av_vel      .~ 0.5 * (vel1 + vel2)
#         norm1        ~ norm(segment)
#         unit_vector .~ segment / norm1
#         k2           ~ 0.1 * k  # compression stiffness tether segments
#         c1           ~ 6.0 * c  # damping kite segments
#         c2           ~ 0.05 * c  # damping perpendicular
#         spring_vel   ~ rel_vel ⋅ unit_vector
#     ]

#     for j in 1:3
#         eqs = [
#             eqs
#             spring_force[j] ~
#                 ((k  * (l_0 - norm1) - c * spring_vel) * unit_vector[j]) * (1 + smooth_sign_ϵ(norm1 - l_0; s.ϵ)) / 2 +
#                 ((k2 * (l_0 - norm1) - c * spring_vel) * unit_vector[j]) * (1 - smooth_sign_ϵ(norm1 - l_0; s.ϵ)) / 2
#         ]
#     end
#     eqs = [
#         eqs
#         v_apparent       ~ v_wind_tether - av_vel
#         i >= s.i_A ?
#             area             ~ norm1 * d_tether * 10 : # 10 is the number of parallel lines in the bridle system
#             area             ~ norm1 * d_tether * (1 + (i%3 == 0)) # double area for middle tether
#         v_app_perp       ~ v_apparent - (v_apparent ⋅ unit_vector) * unit_vector
#         half_drag_force .~ (0.25 * rho * s.set.cd_tether * norm(v_app_perp) * area) .* v_app_perp
#     ]

#     for j in 1:3
#         force_eqs[j, p1] = 
#             (force[j, p1] ~ force_eqs[j, p1].rhs + (half_drag_force[j] + spring_force[j]))
#         force_eqs[j, p2] = 
#             (force[j, p2] ~ force_eqs[j, p2].rhs + (half_drag_force[j] - spring_force[j]))
#     end
    
#     return eqs, force_eqs
# end


function deform_bridle(bridle_pos, static_pos, angle)
end

abstract type AbstractPoint end

struct Point <: AbstractPoint
    idx::Int
    pos::Union{Vector{SimFloat}, Nothing} # pos relative to kite COM in body frame
end

struct WinchPoint <: AbstractPoint
    idx::Int
    pos::Union{Vector{SimFloat}, Nothing} # pos relative to kite COM in body frame
    winch::AbstractWinchModel
end

struct KitePoint <: AbstractPoint
    idx::Int
    pos::Union{Vector{SimFloat}, Nothing} # pos relative to kite COM in body frame
    fixed_pos::Vector{SimFloat} # position in body frame which the point rotates around under kite deformation
    local_y::MVec3 # y-axis in body frame which the point rotates around under kite deformation
end

struct Segment
    idx::Int
    points::Tuple{Int, Int}
    l0::Union{SimFloat, Nothing}
    type::String
end

struct Pulley
    idx
    segments::Tuple{Int, Int}
    sum_length::Union{SimFloat, Nothing}
end

struct PointMassSystem
    points::Vector{AbstractPoint}
    segments::Vector{Segment}
    pulleys::Vector{Pulley}
end

function create_sys!(s::KPSQ; init=false)

    points = Point[]
    segments = Segment[]
    pulleys = Pulley[]
    
    point_idx = 1
    tether_idx = 1
    pulley_idx = 1
    
    bridle_gammas = find_bridle_gammas!(s, zeros(4))

    function create_bridle(gammas)
        i_pnt = length(points) # last point idx
        i_seg = length(segments) # last segment idx
        i_pul = length(pulleys) # last pulley idx

        i = 1
        for gamma in gammas # 2 gammas
            le_pos = [le_interp[i](gamma) for i in 1:3]
            chord = [te_interp[i](gamma) for i in 1:3] .- le_pos
            for frac in bridle_fracs # 4 fracs
                pos = le_pos .+ chord .* frac
                points = [points; KitePoint(i+i_pnt, pos, )]
                i += 1
            end
        end

        points = [
            points
            Point(9+i_pnt, "bridle", [0, 0, 0])
            Point(10+i_pnt, "bridle", [0, 0, 0])
            Point(11+i_pnt, "bridle", [0, 0, 0])
            Point(12+i_pnt, "bridle", [0, 0, 0])

            Point(13+i_pnt, "bridle", [0, 0, -1])

            Point(14+i_pnt, "bridle", [0, 0, -2])
            Point(15+i_pnt, "bridle", [0, 0, -2])

            Point(16+i_pnt, "bridle", [0, 0, -5])
            Point(17+i_pnt, "bridle", [0, 0, -5])
        ]
        segments = [
            segments
            Segment(1+i_seg, (1+i_pnt, 9+i_pnt), 2, "bridle")
            Segment(2+i_seg, (2+i_pnt, 10+i_pnt), 2, "bridle")
            Segment(3+i_seg, (3+i_pnt, 11+i_pnt), 2, "bridle")
            Segment(4+i_seg, (4+i_pnt, 12+i_pnt), 2, "bridle")

            Segment(5+i_seg, (5+i_pnt, 9+i_pnt), 2, "bridle")
            Segment(6+i_seg, (6+i_pnt, 10+i_pnt), 2, "bridle")
            Segment(7+i_seg, (7+i_pnt, 11+i_pnt), 2, "bridle")
            Segment(8+i_seg, (8+i_pnt, 12+i_pnt), 2, "bridle")

            Segment(9+i_seg, (9+i_pnt, 14+i_pnt), 2, "bridle")
            Segment(10+i_seg, (10+i_pnt, 13+i_pnt), 1, "bridle")
            Segment(11+i_seg, (11+i_pnt, 15+i_pnt), 2, "bridle")
            Segment(12+i_seg, (12+i_pnt, 16+i_pnt), 5, "bridle")
            
            Segment(13+i_seg, (13+i_pnt, 14+i_pnt), 1, "bridle")
            Segment(14+i_seg, (13+i_pnt, 15+i_pnt), 1, "bridle")
            
            Segment(15+i_seg, (14+i_pnt, 16+i_pnt), 3, "bridle")
            Segment(16+i_seg, (15+i_pnt, 16+i_pnt), 3, "bridle")
            Segment(17+i_seg, (15+i_pnt, 17+i_pnt), 3, "bridle")
        ]
        pulleys = [
            pulleys
            Pulley(1+i_pul, (13+i_seg, 14+i_seg), nothing)
            Pulley(2+i_pul, (16+i_seg, 17+i_seg), nothing)
        ]
        return 16+i_pnt, 17+i_pnt
    end

    function create_tether(attach_idx)
        l0 = s.set.l_tether / s.set.segments
        for i in 1:s.set.segments
            i_pnt = length(points) # last point idx
            i_seg = length(segments) # last segment idx
            if i == s.set.segments
                points = [points; WinchPoint(1+i_pnt, [0, 0, -5 - i*l0])]
                segments = [segments; Segment(1+i_seg, (i_pnt, 1+i_pnt), l0, "tether")]
            else
                points = [points; Point(1+i_pnt, "tether", [0, 0, -5 - i*l0])]
                segments = [segments; Segment(1+i_seg, (i_pnt, 1+i_pnt), l0, "tether")]
            end
        end
    end

    bridle_attach_idxs = zeros(Int16, 4)
    bridle_attach_idxs[1:2] .= create_bridle(bridle_gammas[1:2])
    bridle_attach_idxs[3:4] .= create_bridle(bridle_gammas[3:4])
    winch_idxs = create_tether.(bridle_attach_idxs)
    @show winch_idxs

    eqs = []
    function force_eqs!()
        pulley_damping = 10
        for pulley in pulleys
            segment1, segment2 = segments[pulley.segments[1]], segments[pulley.segments[2]]
            pulley.sum_length = segment1.l0 + segment2.l0
            mass = pulley.sum_length * segment1.mass_per_meter
            eqs = [
                eqs
                pulley_force[pulley.idx]    ~ spring_force[pulley.segments[1]] - spring_force[pulley.segments[2]]
                pulley_acc[pulley.idx]      ~ pulley_force[pulley.idx] / mass - pulley_damping * pulley_vel[pulley.idx]
            ]
        end

        @variables begin
            segment(t)[1:3, eachindex(segments)]
            unit_vector(t)[1:3, eachindex(segments)]
            l_spring(t), c_spring(t), damping(t), m_tether_particle(t)
            len(t)[eachindex(segments)]
            l0(t)[eachindex(segments)]
            rel_vel(t)[1:3, eachindex(segments)]
            spring_vel(t)[eachindex(segments)]
            spring_force(t)[eachindex(segments)]
            spring_force_vec(t)[1:3, eachindex(segments)]
        end
        for segment in segments
            found = false
            if segment.type == "bridle"
                for pulley in pulleys
                    if segment.idx == pulley.segments[1] # each segment should only be part of one pulley
                        eqs = [
                            eqs
                            l0[segment.idx] ~ pulley_l0[pulley.idx]
                        ]
                        found = true
                        break
                    elseif segment.idx == pulley.segments[2]
                        eqs = [
                            eqs
                            l0[segment.idx] ~ pulley.sum_length - pulley_l0[pulley.idx]
                        ]
                        found = true
                        break
                    end
                end
            end
            if !found
                eqs = [
                    eqs
                    l0[segment.idx] ~ segment.l0
                ]
            end

            p1, p2 = segment.points[1], segment.points[2]

            (segment.type == "bridle") && diameter = s.bridle_tether_diameter
            (segment.type == "power") && diameter = s.power_tether_diameter
            (segment.type == "steering") && diameter = s.steering_tether_diameter

            stiffness = s.set.e_tether * (diameter/2000)^2 * pi
            (segment.type == "bridle") && compression_frac = 1.0
            (segment.type == "power") && compression_frac = 0.1
            (segment.type == "steering") && compression_frac = 0.1
            
            damping = (s.set.damping / s.set.c_spring) * stiffness
            @show damping
            (segment.type == "bridle") && damping = 10damping

            eqs = [
                eqs

                # spring force equations
                segment[:, segment.idx]      ~ pos[:, p2] - pos[:, p1]
                len[segment.idx]             ~ norm(segment[:, segment.idx])
                unit_vector[:, segment.idx]  ~ segment[:, segment.idx]/len[segment.idx]
                rel_vel[:, segment.idx]      ~ vel[:, p1] - vel[:, p2]
                spring_vel[segment.idx]      ~ rel_vel[:, segment.idx] ⋅ unit_vector[:, segment.idx]
                spring_force[segment.idx]    ~ (stiffness * segment.l0 * (len[segment.idx] - l0[segment.idx]) - 
                                            damping * segment.l0 * spring_vel[segment.idx])
                stiffness[segment.idx]       ~ ifelse(len[segment.idx] > l0[segment.idx],
                                            stiffness / len[segment.idx],
                                            compression_frac * stiffness / len[segment.idx]
                )
                damping[segment.idx]         ~ damping / len[segment.idx]
                spring_force ~  (stiffness[segment.idx] * segment.l0 * (len[segment.idx] - l0[segment.idx]) - 
                                damping[segment.idx] * segment.l0 * spring_vel[segment.idx])
                spring_force_vec[:, segment.idx]  ~ spring_force[segment.idx] * unit_vector[:, segment.idx]

                # drag force equations
                height[segment.idx]          ~ max(0.0, 0.5(pos[:, p1][3] + pos[:, p2][3]))
                segment_vel[:, segment.idx]   ~ 0.5(vel[:, p1] + vel[:, p2])
                segment_rho[segment.idx]      ~ calc_rho(s.am, height[i])
                wind_vel[:, segment.idx]     ~ AtmosphericModels.calc_wind_factor(s.am, height[i], s.set.profile_law) * wind_vec_gnd
                app_wind_vel[:, segment.idx] ~ wind_vel - segment_vel
                area[segment.idx]            ~ len[segment.idx] * segment.diameter
                app_perp_vel[:, segment.idx] ~ app_wind_vel[:, segment.idx] - 
                                            (app_wind_vel[:, segment.idx] ⋅ unit_vector[:, segment.idx]) * unit_vector[:, segment.idx]
                drag_force                  ~ (0.5 * tether_rho[segment.idx] * s.set.cd_tether * norm(app_wind_vel[:, segment.idx]) * 
                                            area[segment.idx]) * app_perp_vel[:, segment.idx]
            ]
        end
    
        for point in points
            if point isa WinchPoint
                eqs = [
                    eqs
                    force[:, point.idx]  ~ zeros(3)
                    pos[:, point.idx]    ~ zeros(3)
                    vel[:, point.idx]    ~ zeros(3)
                    acc[:, point.idx]    ~ zeros(3)
                ]
            elseif point isa KitePoint
                eqs = [
                    eqs
                    force[:, point.idx]  ~ zeros(3)
                    pos[:, point.idx]    ~ deform_bridle(point.pos, static_pos, alpha[1])
                    vel[:, point.idx]    ~ zeros(3)
                    acc[:, point.idx]    ~ zeros(3)
                ]
            else
                # segment - inverted
                F::Vector{Num} = zeros(Num, 3)
                mass_per_meter = s.set.rho_tether * π * (segment.diameter/2000)^2    
                mass = 0.
                for segment in segments
                    if point.idx in segment.points
                        inverted = segment.points[2] == point.idx
                        if inverted
                            F .-= spring_force_vec[:, segment.idx]
                        else
                            F .+= spring_force_vec[:, segment.idx]
                        end
                        mass += mass_per_meter * segment.l0 / 2
                        F .+= 0.5drag_force[:, segment.idx]
                    end
                end
                eqs = [
                    eqs
                    force[:, point.idx]  ~ F
                    D(pos[:, point.idx]) ~ vel[:, point.idx]
                    D(vel[:, point.idx]) ~ acc[:, point.idx]
                    acc[:, point.idx]    ~ force[:, point.idx] / mass + G_EARTH
                ]
            end
        end
    end


    if s.torque_control
        [s.motors[i] = TorqueControlledMachine(s.set) for i in 1:3]
    else
        [s.motors[i] = AsyncMachine(s.set) for i in 1:3]
    end
    @parameters begin
        measured_wind_dir_gnd = s.measure.wind_dir_gnd
        measured_sphere_pos[1:2, 1:2] = s.measure.sphere_pos
        measured_sphere_vel[1:2, 1:2] = s.measure.sphere_vel
        measured_sphere_acc[1:2, 1:2] = s.measure.sphere_acc
        measured_tether_length[1:3] = s.measure.tether_length
        measured_tether_vel[1:3]    = s.measure.tether_vel
        measured_tether_acc[1:3]    = s.measure.tether_acc
    end
    @variables begin
        set_values(t)[1:3] # left right middle

        # potential differential variables
        pos(t)[1:3, 1:s.i_C] # xyz pos of left right middle tether
        vel(t)[1:3, 1:s.i_C] 
        acc(t)[1:3, 1:s.i_A-1]
        kite_pos(t)[1:3] # xyz pos of kite in world frame
        kite_vel(t)[1:3]
        kite_acc(t)[1:3]
        distance(t)
        distance_vel(t)
        distance_acc(t)
        wind_scale_gnd(t)
        ω_p(t)[1:3] # turn rate in principal frame
        ω_b(t)[1:3] # turn rate in body frame
        α_p(t)[1:3] # angular acceleration in principal frame
        α_b(t)[1:3]
        trailing_edge_angle(t)[1:2] # angle left right / C D
        trailing_edge_ω(t)[1:2] # angular vel
        trailing_edge_α(t)[1:2] # angular acc
        tether_length(t)[1:3]
        tether_vel(t)[1:3]
        tether_acc(t)[1:3]

        # rotations and frames
        R_b_w(t)[1:3, 1:3] # rotation of the kite body frame relative to the world frame
        R_p_w(t)[1:3, 1:3] # rotation of the kite principal frame relative to the world frame
        Q_p_w(t)[1:4] # quaternion orientation of the kite principal frame relative to the world frame
        Q_b_w(t)[1:4] # quaternion orientation of the kite body frame relative to the world frame
        Q_vel(t)[1:4] # quaternion rate of change
        e_x(t)[1:3]
        e_y(t)[1:3]
        e_z(t)[1:3]
        e_r_D(t)[1:3]
        e_r_E(t)[1:3]
        e_te_A(t)[1:3]
        e_te_B(t)[1:3]

        # rest: forces, torques, vectors and scalar values
        torque_p(t)[1:3] # torque in principal frame
        torque_b(t)[1:3] # torque in body frame
        segment_length(t)[1:3]
        mass_tether_particle(t)[1:3]
        damping(t)[1:3]
        c_spring(t)[1:3]
        total_kite_force(t)[1:3]
        tether_force(t)[1:3] # normalized tether forces at the winch
        force(t)[1:3, 1:s.i_C]
        rho_kite(t)
        norm1(t)[1:s.i_C-3]
        wind_vec_gnd(t)[1:3]
        v_wind_kite(t)[1:3]
    end
    # Collect the arrays into variables
    pos = collect(pos)
    vel = collect(vel)
    acc = collect(acc)

    eqs = []
    force_eqs = SizedArray{Tuple{3, s.i_C}, Symbolics.Equation}(undef)
    force_eqs[:, :] .= (force[:, :] .~ 0)

    te_length = s.kite_length_D/4

    function diff_eqs!()
        Q_b_p = quaternion_conjugate(s.Q_p_b)
        if !init
            Ω = [0       -ω_p[1]  -ω_p[2]  -ω_p[3];
                ω_p[1]   0        ω_p[3]  -ω_p[2];
                ω_p[2]  -ω_p[3]   0        ω_p[1];
                ω_p[3]   ω_p[2]  -ω_p[1]   0]
            eqs = [
                eqs
                [D(Q_p_w[i]) ~ Q_vel[i] for i in 1:4]
                [Q_vel[i] ~ 0.5 * sum(Ω[i, j] * Q_p_w[j] for j in 1:4) for i in 1:4]
                [R_b_w[:, i] ~ quaternion_to_rotation_matrix(quaternion_multiply(Q_p_w, Q_b_p))[:, i] for i in 1:3] # https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Performance_comparisons
                [R_p_w[:, i] ~ quaternion_to_rotation_matrix(Q_p_w)[:, i] for i in 1:3]
                D(ω_p) ~ α_p
                ω_b ~ s.R_b_p' * ω_p
                α_p[1] ~ (torque_p[1] + (s.I_p[2] - s.I_p[3]) * ω_p[2] * ω_p[3]) / s.I_p[1]
                α_p[2] ~ (torque_p[2] + (s.I_p[3] - s.I_p[1]) * ω_p[3] * ω_p[1]) / s.I_p[2]
                α_p[3] ~ (torque_p[3] + (s.I_p[1] - s.I_p[2]) * ω_p[1] * ω_p[2]) / s.I_p[3]
                α_b ~ s.R_b_p' * α_p
                
                [D(kite_pos[i]) ~ kite_vel[i] for i in 1:3]
                [D(kite_vel[i]) ~ kite_acc[i] for i in 1:3]
                kite_acc        ~ total_kite_force / s.set.mass
                distance            ~ norm(kite_pos - pos[:, 3])
                distance_vel        ~ kite_vel ⋅ normalize(kite_pos - pos[:, 3])
                distance_acc        ~ kite_acc ⋅ normalize(kite_pos - pos[:, 3])    
                D(wind_scale_gnd) ~ 0

                [pos[:, i]              .~ 0.0 for i in 1:3]
                [D.(pos[:, i])          .~ vel[:, i] for i in 4:s.i_A-1]
                D(trailing_edge_angle)   ~ trailing_edge_ω
                [vel[:, i]              .~ 0.0 for i in 1:3]
                [D.(vel[:, i])          .~ acc[:, i] for i in 4:s.i_A-1]
                D(trailing_edge_ω)       ~ trailing_edge_α
                D.(tether_length)       .~ tether_vel
                D.(tether_vel)          .~ tether_acc
            ]
        else
            idamp = 10
            scale = 1e-2

            # no movement around body z axis
            Ω = [0       -ω_b[1]  -ω_b[2]  -0     ;
                ω_b[1]   0        0       -ω_b[2];
                ω_b[2]  -0        0        ω_b[1];
                0        ω_b[2]  -ω_b[1]   0     ]
            
            # from measurements
            elevation     = mean(measured_sphere_pos[1, :])
            azimuth       = mean(measured_sphere_pos[2, :])
            elevation_vel = mean(measured_sphere_vel[1, :])
            azimuth_vel   = mean(measured_sphere_vel[2, :])
            elevation_acc = mean(measured_sphere_acc[1, :])
            azimuth_acc   = mean(measured_sphere_acc[2, :])

            r = (measured_sphere_pos[:, 2] - measured_sphere_pos[:, 1]) / 2
            perp_r = [-r[2], r[1]]
            rot_vel = (measured_sphere_vel[:, 1] - measured_sphere_vel[:, 2]) ⋅ (perp_r / norm(r))
            rot_acc = (measured_sphere_acc[:, 1] - measured_sphere_acc[:, 2]) ⋅ (perp_r / norm(r))
            ω_z = rot_vel / norm(r)
            α_z = rot_acc / norm(r)

            angular_acc = measured_tether_acc / s.set.drum_radius
            net_torque = angular_acc * s.set.inertia_total
            measured_tether_force = (net_torque - set_values) / s.set.drum_radius

            # scaled down, partially fixed, and damped equations to find unmeasured variables
            eqs = [
                eqs
                [D(Q_b_w[i]) ~ Q_vel[i] for i in 1:4]
                Q_p_w ~ quaternion_multiply(Q_b_w, quaternion_conjugate(Q_b_p))
                [Q_vel[i] ~ 0.5 * sum(Ω[i, j] * Q_b_w[j] for j in 1:4) for i in 1:4]
                [R_b_w[:, i] ~ quaternion_to_rotation_matrix(Q_b_w)[:, i] for i in 1:3] # https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Performance_comparisons
                [R_p_w[:, i] ~ quaternion_to_rotation_matrix(Q_p_w)[:, i] for i in 1:3]
                D(ω_b[1:2]) ~ scale * α_b[1:2] - idamp * ω_b[1:2]
                ω_b[3] ~ ω_z
                ω_p ~ s.R_b_p * ω_b
                torque_b ~ s.R_b_p' * torque_p
                α_b[1] ~ (torque_b[1]) / s.I_b[1]
                α_b[2] ~ (torque_b[2]) / s.I_b[2]
                α_b[3] ~ α_z

                kite_pos        ~ rotate_around_z(rotate_around_y([distance, 0, 0], -elevation), azimuth)
                kite_vel        ~ distance_vel * normalize(kite_pos - pos[:, 3]) + 
                                    rotate_around_z(rotate_around_y([0, azimuth_vel * distance, elevation_vel * distance], -elevation), azimuth)
                kite_acc        ~ distance_acc * normalize(kite_pos - pos[:, 3]) + 
                                    rotate_around_z(rotate_around_y([0, azimuth_acc * distance, elevation_acc * distance], -elevation), azimuth)
                D(distance)     ~ distance_vel
                D(distance_vel) ~ distance_acc - idamp * distance_vel
                distance_acc    ~ scale * (measured_tether_acc[3] - tether_acc[3])
                D(wind_scale_gnd)  ~ 100scale * (mean(measured_tether_force[1:2]) - mean(tether_force[1:2]))

                [pos[:, i]              .~ 0.0 for i in 1:3]
                [D.(pos[:, i])          .~ vel[:, i] for i in 4:s.i_A-1]
                D(trailing_edge_angle)   ~ trailing_edge_ω
                [vel[:, i]              .~ 0.0 for i in 1:3]
                [D.(vel[:, i])          .~ scale * acc[:, i] - idamp * vel[:, i] for i in 4:s.i_A-1] # TODO: ADD CENTRIFUGAL FORCE DUE TO KITE ROTATION
                D(trailing_edge_ω)       ~ 0.1scale * trailing_edge_α - idamp * trailing_edge_ω
                tether_length           ~ measured_tether_length
                tether_vel              ~ measured_tether_vel
            ]
        end
        return nothing
    end
    function scalar_eqs!()
        # last tether point vel without rotational trailing edge vel
        kite_vel_A = kite_vel + R_b_w * (ω_b × s.pos_A_b)
        kite_vel_B = kite_vel + R_b_w * (ω_b × s.pos_B_b)
        kite_vel_C = kite_vel + R_b_w * (ω_b × s.pos_C_b)

        mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2

        eqs = [
            eqs

            # last tether point pos
            pos[:, s.i_C]    ~ kite_pos + R_b_w * s.pos_C_b
            pos[:, s.i_A]    ~ pos[:, s.i_C] + e_y * s.pos_D_b[2] + e_x * te_length * cos(trailing_edge_angle[1]) + # TODO: fix pos
                e_r_D * te_length * sin(trailing_edge_angle[1])
            pos[:, s.i_B]    ~ pos[:, s.i_C] + e_y * s.pos_E_b[2] + e_x * te_length * cos(trailing_edge_angle[2]) + 
                e_r_E * te_length * sin(trailing_edge_angle[2])

            # last tether point vel with rotational trailing edge vel
            vel[:, s.i_C]   ~ kite_vel_C
            vel[:, s.i_A]   ~ kite_vel_A + e_x * te_length * cos(trailing_edge_ω[1]) + e_r_D * te_length * sin(trailing_edge_ω[1])
            vel[:, s.i_B]   ~ kite_vel_B + e_x * te_length * cos(trailing_edge_ω[2]) + e_r_E * te_length * sin(trailing_edge_ω[2])

            segment_length          ~ tether_length  ./ s.set.segments
            mass_tether_particle    ~ mass_per_meter .* segment_length
            # damping                 ~ [s.damping / segment_length[1], s.damping / segment_length[2], s.damping*2 / segment_length[3]]
            damping                 ~ s.damping ./ segment_length
            # c_spring                ~ [s.c_spring / segment_length[1], s.c_spring / segment_length[2], s.c_spring*2 / segment_length[3]]
            c_spring                ~ s.c_spring ./ segment_length
            e_x     ~ R_b_w * [1, 0, 0]
            e_y     ~ R_b_w * [0, 1, 0]
            e_z     ~ R_b_w * [0, 0, 1]
            e_r_D   ~ rotate_v_around_k(e_z, e_x, 0.5π + s.γ_D)
            e_r_E   ~ rotate_v_around_k(e_z, e_x, 1.5π - s.γ_D)
            e_te_A  ~ -e_x * sin(trailing_edge_angle[1]) + e_r_D * cos(trailing_edge_angle[1])
            e_te_B  ~ -e_x * sin(trailing_edge_angle[2]) + e_r_E * cos(trailing_edge_angle[2])
            rho_kite        ~ calc_rho(s.am, pos[3,s.i_A])
            wind_vec_gnd ~ wind_scale_gnd * rotate_around_z([1, 0, 0], measured_wind_dir_gnd)
            v_wind_kite     ~ AtmosphericModels.calc_wind_factor(s.am, kite_pos[3], s.set.profile_law) * wind_vec_gnd
        ]

        @variables begin
            heading_y(t)
            power_angle(t) # average flap angle
            power_vel(t)
            steering_angle(t) # difference between left and right flap angle
            steering_vel(t)
            tether_diff(t)
            tether_diff_vel(t)
            set_diff(t)
            azimuth(t)
            azimuth_vel(t)
            azimuth_acc(t)
            elevation(t)
            elevation_vel(t)
            elevation_acc(t)
            x_acc(t)
            y_acc(t)
            left_diff(t)
            right_diff(t)
            sphere_pos(t)[1:2, 1:2] # TODO: add equations for these, and think of a good measurement system, maybe rolling window?
            sphere_vel(t)[1:2, 1:2]
            sphere_acc(t)[1:2, 1:2]
        end

        x, y, z = kite_pos
        x´, y´, z´ = kite_vel
        x´´, y´´, z´´ = kite_acc
        eqs = [
            eqs
            tether_force     ~ [norm(force[:, i]) for i in 1:3]
            heading_y       ~ calc_heading_y(-e_x)
            power_angle         ~ (trailing_edge_angle[1] + trailing_edge_angle[2]) / 2
            power_vel           ~ (trailing_edge_ω[1] + trailing_edge_ω[2]) / 2
            steering_angle      ~ trailing_edge_angle[2] - trailing_edge_angle[1]
            steering_vel        ~ trailing_edge_ω[2] - trailing_edge_ω[1]
            tether_diff         ~ tether_length[2] - tether_length[1]
            tether_diff_vel     ~ tether_vel[2] - tether_vel[1]
            set_diff            ~ set_values[2] - set_values[1]

            elevation           ~ atan(z / x)
            # elevation_vel = d/dt(atan(z/x)) = (x*ż' - z*ẋ')/(x^2 + z^2) according to wolframalpha
            elevation_vel       ~ (x*z´ - z*x´) / 
                                    (x^2 + z^2)
            elevation_acc       ~ ((x^2 + z^2)*(x*z´´ - z*x´´) + 2(z*x´ - x*z´)*(x*x´ + z*z´))/(x^2 + z^2)^2
            azimuth             ~ atan(y / x)
            # azimuth_vel = d/dt(atan(y/x)) = (-y*x´ + x*y´)/(x^2 + y^2) # TODO: check if correct
            azimuth_vel         ~ (-y*x´ + x*y´) / 
                                    (x^2 + y^2)
            azimuth_acc         ~ ((x^2 + y^2)*(-y*x´´ + x*y´´) + 2(y*x´ - x*y´)*(x*x´ + y*y´))/(x^2 + y^2)^2
            x_acc               ~ kite_acc ⋅ e_x
            y_acc               ~ kite_acc ⋅ e_y
            left_diff           ~ tether_length[1] - tether_length[3]
            right_diff          ~ tether_length[2] - tether_length[3]
        ]
        return nothing
    end
    # function tether_forces!()
    #     @variables begin
    #         height(t)[1:s.i_C-3]
    #         rho(t)[1:s.i_C-3]
    #         v_wind_tether(t)[1:3, 1:s.i_C-3]
    #         l_0(t)[1:s.i_C-3]
    #         k(t)[1:s.i_C-3]
    #         c(t)[1:s.i_C-3]
    #         segment(t)[1:3, 1:s.i_C-3]
    #         rel_vel(t)[1:3, 1:s.i_C-3]
    #         av_vel(t)[1:3, 1:s.i_C-3] 
    #         unit_vector(t)[1:3, 1:s.i_C-3]
    #         k2(t)[1:s.i_C-3]
    #         c1(t)[1:s.i_C-3]
    #         c2(t)[1:s.i_C-3]
    #         spring_vel(t)[1:s.i_C-3]
    #         spring_force(t)[1:3, 1:s.i_C-3]
    #         v_apparent(t)[1:3, 1:s.i_C-3]
    #         area(t)[1:s.i_C-3]
    #         v_app_perp(t)[1:3, 1:s.i_C-3]
    #         half_drag_force(t)[1:3, 1:s.i_C-3]
    #     end
        
    #     for i in 1:s.i_C-3
    #         p1 = i  # First point nr.
    #         p2 = i+3
    #         eqs = [
    #             eqs
    #             height[i]           ~ max(0.0, 0.5 * (pos[:, p1][3] + pos[:, p2][3]))
    #             rho[i]              ~ calc_rho(s.am, height[i])
    #             v_wind_tether[:, i] ~ AtmosphericModels.calc_wind_factor(s.am, height[i], s.set.profile_law) * wind_vec_gnd
    #         ]
    
    #         eqs, force_eqs = calc_particle_forces!(s, eqs, force_eqs, force, p1, p2, pos[:, p1], pos[:, p2], vel[:, p1], 
    #                           vel[:, p2], segment_length, c_spring, damping, rho[i], i, l_0[i], k[i], c[i], segment[:, i], 
    #                           rel_vel[:, i], av_vel[:, i], norm1[i], unit_vector[:, i], k2[i], c1[i], c2[i], spring_vel[i],
    #                           spring_force[:, i], v_apparent[:, i], v_wind_tether[:, i], area[i], v_app_perp[:, i],
    #                           half_drag_force[:, i])
    #     end
    #     return nothing
    # end
    function kite_forces!()
        n = s.set.aero_surfaces

        # integrating loop variables, iterating over 2n segments
        @variables begin
            aero_pos(t)[1:3, 1:2n]
            e_r(t)[1:3, 1:2n]
            e_te(t)[1:3, 1:2n] # clockwise trailing edge of flap vector
            y_l(t)[1:2n]
            v_a(t)[1:3, 1:2n]
            e_drift(t)[1:3, 1:2n]
            v_a_xr(t)[1:3, 1:2n]
            aoa(t)[1:n*2]
            seg_cl(t)[1:n*2]
            seg_cd(t)[1:n*2]
            seg_L(t)[1:3, 1:2n]
            seg_D(t)[1:3, 1:2n]
            seg_aero_force(t)[1:3, 1:2n]
            seg_g_force(t)[1:3, 1:2n]
            seg_te_force(t)[1:3, 1:2n]
            seg_trailing_edge_angle(t)[1:2n] # flap angle relative to -e_x, e_z
            ram_force(t)[1:2n]
            te_force(t)[1:2n]
            seg_aero_torque_p(t)[1:3, 1:2n]
            seg_gravity_torque_p(t)[1:3, 1:2n]
            C_torque_p(t)[1:3] # torque on point C in principal reference frame
            seg_vel(t)[1:3, 1:2n]
            F(t)[1:3, 1:2n]
            r(t)[1:3, 1:2n]
            L_C(t)[1:3]
            L_D(t)[1:3]
            D_C(t)[1:3]
            D_D(t)[1:3]
            F_te_C(t)[1:3] # trailing edge clockwise force
            F_te_D(t)[1:3]
        end

        s.γ_l         = π/2 - s.set.width/2/s.set.radius
        γ_middle    = π/2
        dγ          = (γ_middle - s.γ_l) / n
        γ_m_l = π/2 - s.set.min_steering_line_distance/s.set.radius # end of steering lines on left side
        γ_m_r = π/2 + s.set.min_steering_line_distance/s.set.radius # end of steering lines on right side
        
        for i in 1:2n
            if i <= n
                γ = s.γ_l + -dγ/2 + i * dγ
                kite_length = s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (γ - s.γ_l) / (π/2 - s.γ_l) # TODO: kite length gets less with flap turning
            else
                γ = pi - (s.γ_l + -dγ/2 + (i-n) * dγ)
                kite_length = s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (π - s.γ_l - γ) / (π/2 - s.γ_l)
            end
            seg_flap_height = kite_length * s.set.flap_height

            eqs = [
                eqs
                e_r[:, i]       ~ rotate_v_around_k(e_z, e_x, 0.5π + γ)
                seg_vel[:, i]   ~ kite_vel + R_b_w * (ω_b × s.seg_cop_pos_b[:, i])
                v_a[:, i]       ~ v_wind_kite .- seg_vel[:, i]
                e_drift[:, i]   ~ (e_x × e_r[:, i])
                v_a_xr[:, i]    ~ v_a[:, i] .- (v_a[:, i] ⋅ e_drift[:, i]) .* e_drift[:, i]

                γ < γ_m_l ?
                    seg_trailing_edge_angle[i] ~ trailing_edge_angle[1]  :
                γ > γ_m_r ?
                    seg_trailing_edge_angle[i] ~ trailing_edge_angle[2] :
                    seg_trailing_edge_angle[i] ~ ((trailing_edge_angle[2] - trailing_edge_angle[1]) / (γ_m_r - γ_m_l) * (γ - γ_m_l) + trailing_edge_angle[1])

                aoa[i]      ~ -asin2(normalize(v_a_xr[:, i]) ⋅ e_r[:, i]) + deg2rad(s.set.alpha_zero)
                seg_cl[i]   ~ sym_interp(s.cl_interp, aoa[i], seg_trailing_edge_angle[i])
                seg_cd[i]   ~ sym_interp(s.cd_interp, aoa[i], seg_trailing_edge_angle[i])

                seg_L[:, i] ~ 0.5 * rho_kite * (norm(v_a_xr[:, i]))^2 * s.set.radius * dγ * kite_length * seg_cl[i] * 
                                    normalize(v_a_xr[:, i] × e_drift[:, i])
                seg_D[:, i] ~ 0.5 * rho_kite * norm(v_a_xr[:, i]) * s.set.radius * dγ * kite_length * seg_cd[i] *
                                    v_a_xr[:, i]

                e_te[:, i] ~ -e_x * sin(seg_trailing_edge_angle[i]) + e_r[:, i] * cos(seg_trailing_edge_angle[i])
                ram_force[i] ~ smooth_sign_ϵ(deg2rad(s.set.alpha_zero) - seg_trailing_edge_angle[i]; s.ϵ) *
                            rho_kite * norm(v_a[:, i])^2 * seg_flap_height * s.set.radius * dγ * (seg_flap_height/2) / (kite_length/4)
                te_force[i] ~ 0.5 * rho_kite * (norm(v_a_xr[:, i]))^2 * s.set.radius * dγ * kite_length * 
                            sym_interp(s.c_te_interp, aoa[i], seg_trailing_edge_angle[i])
                seg_te_force[:, i] ~ (ram_force[i] + te_force[i]) * e_te[:, i]
                
                seg_aero_force[:, i] ~ seg_L[:, i] + seg_D[:, i]
                seg_g_force[:, i] ~ [0.0, 0.0, -G_EARTH * s.seg_mass[i]]
                seg_aero_torque_p[:, i] ~ s.seg_cop_pos_p[:, i] × (R_p_w' * seg_aero_force[:, i])
                seg_gravity_torque_p[:, i] ~ s.seg_com_pos_p[:, i] × (R_p_w' * seg_g_force[:, i])
            ]
        end

        eqs = [
            eqs
            total_kite_force ~ [sum(seg_aero_force[i, :]) + sum(seg_g_force[i, :]) + force[i, s.i_C] for i in 1:3]
            C_torque_p ~ s.pos_C_p × (R_p_w' * force[:, s.i_C])
            torque_p ~ [sum(seg_aero_torque_p[i, :]) + sum(seg_gravity_torque_p[i, :]) + C_torque_p[i] for i in 1:3]
            F_te_C ~ [sum(seg_te_force[i, 1:n]) for i in 1:3]
            F_te_D ~ [sum(seg_te_force[i, n+1:2n]) for i in 1:3]
        ]

        force_eqs[:,s.i_A] .= (force[:, s.i_A] .~ F_te_C + ([0.0, 0.0, -G_EARTH*(s.set.mass/8)]) * (s.set.mass/8))
        force_eqs[:,s.i_B] .= (force[:, s.i_B] .~ F_te_D + ([0.0, 0.0, -G_EARTH*(s.set.mass/8)]) * (s.set.mass/8))
        return nothing
    end

    diff_eqs!()
    scalar_eqs!()
    tether_forces!()
    kite_forces!()
    
    if s.torque_control
        eqs = vcat(eqs, tether_acc .~ [calc_acc_torque(s.motors[i], tether_vel[i], norm(force[:, (i-1)%3+1]),
            set_values[i]) for i in 1:3])
    else
        eqs = vcat(eqs, tether_acc .~ [calc_acc_speed(s.motors[i], tether_vel[i], norm(force[:, (i-1)%3+1]),
            set_values[i]) for i in 1:3])
    end
    for i in 1:3
        eqs = vcat(eqs, vcat(force_eqs[:, i]))
        eqs = vcat(eqs, acc[:, i] .~ 0)
    end
    for i in 4:s.i_A-1
        eqs = vcat(eqs, vcat(force_eqs[:, i]))
        eqs = vcat(eqs, acc[:, i] .~ [0.0, 0.0, -G_EARTH] + (force[:, i] / mass_tether_particle[(i-1)%3+1]))
    end

    te_I = (1/3 * (s.set.mass/8) * te_length^2)
    # -damping / I * ω = α_damping
    # solve for c: (c * (k*m/s^2) / (k*m^2)) * (m/s)=m/s^2 in wolframalpha
    # damping should be N*m*s
    rot_damping = 0.1s.damping * te_length

    eqs = [
        eqs
        vcat(force_eqs[:, s.i_A])
        vcat(force_eqs[:, s.i_B])
        vcat(force_eqs[:, s.i_C])
        trailing_edge_α[1] ~ (force[:, s.i_A]) ⋅ e_te_A * te_length / te_I - (rot_damping[1] / te_I) * trailing_edge_ω[1]
        trailing_edge_α[2] ~ (force[:, s.i_B]) ⋅ e_te_B * te_length / te_I - (rot_damping[2] / te_I) * trailing_edge_ω[2]
    ]
    
    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))

    if !init
        discrete_events = [true => [Q_p_w[i] ~ normalize(Q_p_w)[i] for i in 1:4]]
    else
        discrete_events = [true => [Q_b_w[i] ~ normalize(Q_b_w)[i] for i in 1:4]]
    end
    
    @named sys = ODESystem(eqs, t; discrete_events)
    return sys, collect(set_values)
end

"""
The distance/vel/acc of the kite cannot be measured directly, but the average acc of the kite distance is equal to the tether acc.
    So while the set_values input and wind is constant, distance_acc = tether_acc[3]. To accurately describe distance_acc when set_values
    or wind suddenly change, the distance_acc is calculated by combining the simulated distance_acc from timestep-1 to timestep and the
    current measured acc.
    distance_acc = 0.99 * sim_distance_acc + 0.01 * measured_distance_acc

Assume distance_acc = tether_acc[3] for convenience
"""
function model!(s::KPSQ; init=false)
    init_pos!(s)
    
    sys, inputs = create_sys!(s; init)
    # structural_simplify(sys, (inputs, []))
    (sys, _) = structural_simplify(sys, (inputs, []); fully_determined=true)

    if !init
        u0map = [
            [sys.Q_p_w[i] => s.Q_p_w[i] for i in 1:4]
            [sys.ω_p[i] => 0 for i in 1:3]

            [sys.kite_pos[j] => s.kite_pos[j] for j in 1:3]
            [sys.kite_vel[i] => 0.0 for i in 1:3]

            [sys.pos[j, i] => s.pos[j, i] for j in 1:3 for i in 4:s.i_A-1]
            [sys.vel[j, i] => 0 for j in 1:3 for i in 4:s.i_A-1]

            [sys.trailing_edge_angle[i] => s.te_angle[i] for i in 1:2]
            [sys.trailing_edge_ω[i] => 0 for i in 1:2]

            [sys.tether_length[i] => s.measure.tether_length[i] for i in 1:3]
            [sys.tether_vel[j] => 0 for j in 1:3]

            sys.wind_scale_gnd => s.set.v_wind
        ]
    else
        Q_b_w = quaternion_multiply(s.Q_p_w, quaternion_conjugate(s.Q_p_b))
        u0map = [
            [sys.Q_b_w[i] => Q_b_w[i] for i in 1:4]
            [sys.ω_b[i] => 0 for i in 1:2]

            sys.distance => norm(s.kite_pos)
            sys.distance_vel => 0.0

            [sys.pos[j, i] => s.pos[j, i] for j in 1:3 for i in 4:s.i_A-1]
            [sys.vel[j, i] => 0 for j in 1:3 for i in 4:s.i_A-1]

            [sys.trailing_edge_angle[i] => s.te_angle[i] for i in 1:2]
            [sys.trailing_edge_ω[i] => 0 for i in 1:2]

            sys.wind_scale_gnd => s.set.v_wind
        ]
    end
    p0map = [sys.set_values[j] => s.measure.set_values[j] for j in 1:3]
    return sys, u0map, p0map
end

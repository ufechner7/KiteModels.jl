module Environment

using KiteModels, StaticArrays, LinearAlgebra, Parameters
import OrdinaryDiffEqCore: ODEIntegrator
export reset, step, render, Env

const StateVec = MVector{11, Float32}

@with_kw mutable struct Env
    kcu::KCU = KCU(se())
    s::KPSQ = KPSQ(kcu)
    max_render_length::Int = 10000
    i::Int = 1
    logger::Logger = Logger(s.i_A, max_render_length)
    integrator::ODEIntegrator = KiteModels.init_sim!(s; prn=false, torque_control=true)
    sys_state::SysState = SysState(s)
    state::StateVec = zeros(StateVec)
    state_d::StateVec = zeros(StateVec)
    state_dd::StateVec = zeros(StateVec)
    last_state::StateVec = zeros(StateVec)
    last_state_d::StateVec = zeros(StateVec)
    wanted_elevation::Float64 = 0.0
    wanted_azimuth::Float64 = 0.0
    wanted_tether_length::Float64 = 0.0
    max_force::Float64 = 0.0
    rotation::Float64 = 0.0
end

function step(e::Env, reel_out_torques; prn=false)
    reel_out_torques = Vector{Float64}(reel_out_torques)

    old_heading = calc_heading(e.s)
    if prn
        KiteModels.next_step!(e.s, e.integrator; set_values=reel_out_torques)
    else
        redirect_stderr(devnull) do
            redirect_stdout(devnull) do
                KiteModels.next_step!(e.s, e.integrator; set_values=reel_out_torques)
            end
        end
    end
    _calc_rotation(e, old_heading, calc_heading(e.s))
    update_sys_state!(e.sys_state, e.s)
    e.i += 1
    return (e.integrator.last_stepfail, _calc_state(e, e.s))
end

function reset(e::Env, name="sim_log", elevation=0.0, azimuth=0.0, tether_length=50.0, force=5000.0, log=false)
    e.wanted_elevation = Float32(elevation)
    e.wanted_azimuth = Float32(azimuth)
    e.wanted_tether_length = Float32(tether_length)
    e.max_force = Float32(force)
    e.rotation = 0.0
    if length(e.logger) > 1
        name = String(name)
        save_log(e.logger, basename(name))
    end
    update_settings()
    e.logger = Logger(e.s.i_A, e.max_render_length)
    e.integrator = KiteModels.init_sim!(e.s; prn=false, torque_control=true)
    e.sys_state = SysState(e.s)
    e.i = 1
    return _calc_state(e, e.s)
end

function render(e::Env)
    if(e.i <= e.max_render_length)
        log!(e.logger, e.sys_state)
    end
end

function _calc_state(e::Env, s::KPSQ)
    e.state .= vcat(
        _calc_speed_reward(e,s),                # length 1
        calc_orient_quat(s),            # length 4
        s.tether_lengths,                    # length 3 # normalize to min and max 
        
        KiteModels.calc_tether_elevation(s),       # length 1
        KiteModels.calc_tether_azimuth(s),         # length 1
        sum(winch_force(s)),            # length 1
    )
    if e.i == 1
        e.state_d .= zeros(StateVec)
    else
        e.state_d .= (e.state .- e.last_state) / e.s.set.sample_freq
    end
    if e.i <= 2
        e.state_dd .= zeros(StateVec)
    else
        e.state_dd .= (e.state_d .- e.last_state_d) / e.s.set.sample_freq
    end
    e.last_state_d .= e.state_d
    e.last_state .= e.state
    return vcat(e.state, e.state_d, e.state_dd)
end

function _calc_reward(e::Env, s::KPSQ)
    if  s.vel_kite ⋅ s.e_x < 0 ||
        (KiteModels.calc_tether_elevation(s) < e.wanted_elevation ||
        !(-2*π < e.rotation < 2*π) ||
        s.tether_lengths[3] > e.wanted_tether_length*1.5 ||
        s.tether_lengths[3] < e.wanted_tether_length*0.95 ||
        sum(winch_force(s)) > e.max_force)
        return 0.0
    end
    force_component = _calc_force_component(e,s)
    reward = clamp(force_component / e.max_force, 0.0, 1.0) # range [-1, 1] clamped to [0, 1] because 0 is physical minimum
    return reward
end

function _calc_speed_reward(e::Env, s::KPSQ)
    speed = s.vel_kite ⋅ s.e_x
    if  (KiteModels.calc_tether_elevation(s) < e.wanted_elevation ||
        !(-2*π < e.rotation < 2*π) ||
        s.tether_lengths[3] > e.wanted_tether_length*1.5 ||
        s.tether_lengths[3] < e.wanted_tether_length*0.95 ||
        sum(winch_force(s)) > e.max_force)
        return 0.0
    end
    # wanted_minus_z = [cos(e.wanted_elevation)*cos(e.wanted_azimuth), cos(e.wanted_elevation)*-sin(e.wanted_azimuth), sin(e.wanted_elevation)]
    # reward = speed * (-s.e_z ⋅ wanted_minus_z)
    reward = speed
    return reward
end

function _calc_direction_reward(e::Env, s::KPSQ)
    if  s.vel_kite ⋅ s.e_x < 0.0 ||
        (KiteModels.calc_tether_elevation(s) < e.wanted_elevation ||
        !(-2*π < e.rotation < 2*π) ||
        s.tether_lengths[3] > e.wanted_tether_length*1.5 ||
        s.tether_lengths[3] < e.wanted_tether_length*0.95 ||
        sum(winch_force(s)) > e.max_force)
        return 0.01
    end
    wanted_direction = [cos(e.wanted_elevation)*cos(e.wanted_azimuth), cos(e.wanted_elevation)*-sin(e.wanted_azimuth), sin(e.wanted_elevation)]
    reward = normalize(s.pos[6]) ⋅ wanted_direction + 1.0
    return reward
end

function _calc_force_component(e::Env, s::KPSQ)
    wanted_force_vector = [cos(e.wanted_elevation)*cos(e.wanted_azimuth), cos(e.wanted_elevation)*-sin(e.wanted_azimuth), sin(e.wanted_elevation)]
    tether_force = sum(s.winch_forces)
    force_component = tether_force ⋅ wanted_force_vector
    return force_component
end

function _calc_rotation(e::Env, old_heading, new_heading)
    d_rotation = new_heading - old_heading
    if d_rotation > 1
        d_rotation -= 2*pi
    elseif d_rotation < -1
        d_rotation += 2*pi
    end
    e.rotation += d_rotation
    return nothing
end

end
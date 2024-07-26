
using KiteModels, StaticArrays, LinearAlgebra, Serialization
export reset, step, render

const Model = KPS4_3L
kcu = nothing;
s = nothing;
dt = 1/se().sample_freq;
max_render_length = 10000;
i = 1;
logger = nothing;
integrator = nothing;
steady_state_history = load_history()
const StateVec = MVector{11, Float32}
state::StateVec = zeros(StateVec)
state_d::StateVec = zeros(StateVec)
state_dd::StateVec = zeros(StateVec)
last_state::StateVec = zeros(StateVec)
last_state_d::StateVec = zeros(StateVec)
wanted_elevation = 0.0
wanted_azimuth = 0.0
wanted_tether_length = 0.0
max_force = 5000.0
rotation = 0.0

function step(reel_out_speeds; prn=false)
    global i, sys_state, s, integrator
    reel_out_speeds = Vector{Float64}(reel_out_speeds)

    old_heading = calc_heading(s)
    if prn
        KiteModels.next_step!(s, integrator, set_speeds=reel_out_speeds, dt=dt)
    else
        redirect_stdout(devnull) do
            KiteModels.next_step!(s, integrator, set_speeds=reel_out_speeds, dt=dt)
        end
    end
    
    _calc_rotation(old_heading, calc_heading(s))
    i += 1
    return _calc_state(s)
end

function reset(name="sim_log", elevation=0.0, azimuth=0.0, tether_length=50.0, force=5000.0)
    global kcu, s, integrator, steady_state_history, i, sys_state, logger, steady_state_history, wanted_elevation, wanted_azimuth, wanted_tether_length, rotation, max_force

    wanted_elevation = Float32(elevation)
    wanted_azimuth = Float32(azimuth)
    wanted_tether_length = Float32(tether_length)
    max_force = Float32(force)

    rotation = 0.0

    if logger !== nothing
        name = String(name)
        path = joinpath(get_data_path(), name) * ".bin"
        serialize(path, logger)
    end
    
    update_settings()
    kcu = KCU(se());
    s = Model(kcu);
    logger = Vector{typeof(s.pos)}()
    integrator = KiteModels.init_sim!(s, stiffness_factor=0.04, prn=false, steady_state_history=steady_state_history)
    i = 1
    return _calc_state(s)
end

function render()
    global sys_state, logger, i
    if(i <= max_render_length)
        push!(logger, deepcopy(s.pos))
    end
end

function _calc_state(s::KPS4_3L)
    global state, state_d, state_dd, last_state_d, last_state
    _calc_reward(s)
    state .= vcat(
        _calc_reward(s),                # length 1
        calc_orient_quat(s),            # length 4
        s.l_tethers,                    # length 3 # normalize to min and max l_tethers
        KiteModels.calc_tether_elevation(s),       # length 1
        KiteModels.calc_tether_azimuth(s),         # length 1
        sum(winch_force(s)),            # length 1
    )
    if i == 1
        state_d .= zeros(StateVec)
    else
        state_d .= (state .- last_state) / dt
    end
    if i <= 2
        state_dd .= zeros(StateVec)
    else
        state_dd .= (state_d .- last_state_d) / dt
    end
    last_state_d .= state_d
    last_state .= state
    return vcat(state, state_d, state_dd)
end

function _calc_reward(s::KPS4_3L)
    global max_force
    if (KiteModels.calc_tether_elevation(s) < wanted_elevation ||
        !(-2*π < rotation < 2*π) ||
        s.l_tethers[1] > wanted_tether_length*1.5 ||
        s.l_tethers[1] < wanted_tether_length*0.95 ||
        sum(winch_force(s)) > max_force)
        return 0.0
    end
    force_component = _calc_force_component(s)
    reward = clamp(force_component / max_force, 0.0, 1.0) # range [-1, 1] clamped to [0, 1] because 0 is physical minimum
    return reward
end

function _calc_force_component(s::KPS4_3L)
    global wanted_elevation, wanted_azimuth
    wanted_force_vector = [cos(wanted_elevation)*cos(wanted_azimuth), cos(wanted_elevation)*-sin(wanted_azimuth), sin(wanted_elevation)]
    tether_force = sum(s.forces[1:3])
    force_component = tether_force ⋅ wanted_force_vector
    return force_component
end

function _calc_rotation(old_heading, new_heading)
    global rotation
    d_rotation = new_heading - old_heading
    if d_rotation > 1
        d_rotation -= 2*pi
    elseif d_rotation < -1
        d_rotation += 2*pi
    end
    rotation += d_rotation
    return nothing
end


function close()
    save_history(steady_state_history)
end

reset()
for i in 1:1000
    println(
        step([0,0,0])
    )
    render()
end
reset()
close()
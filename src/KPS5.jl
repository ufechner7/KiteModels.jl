#= MIT License

Copyright (c) 2020, 2021, 2022, 2024 Uwe Fechner

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

This model implements a 3D mass-spring system with reel-out. It uses six tether segments (the number can be
configured in the file data/settings.yaml). The kite is modelled using 4 point masses and 3 aerodynamic 
surfaces. The spring constant and the damping decrease with the segment length. The aerodynamic kite forces
are acting on three of the four kite point masses. 

Four point kite model, included from KiteModels.jl.

Scientific background: http://arxiv.org/abs/1406.6218 =#
KiteModels.stiffnea=1000
"""
    mutable struct KPS5{S, T, P, Q, SP} <: AbstractKiteModel

State of the kite power system, using a 4 point kite model. Parameters:
- S: Scalar type, e.g. SimFloat
  In the documentation mentioned as Any, but when used in this module it is always SimFloat and not Any.
- T: Vector type, e.g. MVector{3, SimFloat}
- P: number of points of the system, segments+1
- Q: number of springs in the system, P-1
- SP: struct type, describing a spring
Normally a user of this package will not have to access any of the members of this type directly,
use the input and output functions instead.

$(TYPEDFIELDS)
"""
@with_kw mutable struct KPS5{S, T} <: AbstractKiteModel
    "Reference to the settings struct"
    set::Settings
    "Reference to the KCU model (Kite Control Unit as implemented in the package KitePodModels"
    kcu::KCU
    "Reference to the atmospheric model as implemented in the package AtmosphericModels"
    am::AtmosphericModel = AtmosphericModel()
    "Reference to winch model as implemented in the package WinchModels"
    wm::AbstractWinchModel
    "Iterations, number of calls to the function residual!"
    iter:: Int64 = 0
    "x vector of kite reference frame"
    x::T =                 zeros(S, 3)
    "y vector of kite reference frame"
    y::T =                 zeros(S, 3)
    "z vector of kite reference frame"
    z::T =                 zeros(S, 3)
    sys::Union{ModelingToolkit.ODESystem, Nothing} = nothing
    t_0::Float64 = 0.0
    #iter::Int64 = 0
    prob::Union{OrdinaryDiffEqCore.ODEProblem, Nothing} = nothing
    integrator::Union{OrdinaryDiffEqCore.ODEIntegrator, Nothing} = nothing
    get_state::Function            = () -> nothing
end
function KPS5(kcu::KCU)
    if kcu.set.winch_model == "AsyncMachine"
        wm = AsyncMachine(kcu.set)
    elseif kcu.set.winch_model == "TorqueControlledMachine"
        wm = TorqueControlledMachine(kcu.set)
    end
    # wm.last_set_speed = kcu.set.v_reel_out
    s = KPS5{SimFloat, KVec3}(set=kcu.set, 
             kcu=kcu, wm=wm)    
    return s
end
# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# deriving a constants for points, total segments, and connections
# ---------
function points(s)
    return 5 + s.set.segments 
end
function total_segments(s)
    return 9 + s.set.segments
end
function getconnections(s)
    conn = [(1,2), (2,3), (3,1), (1,4), (2,4), (3,4), (1,5), (2,5), (3,5)]      # connections between KITE points 
    conn = vcat(conn, [(6+i, 6+i+1) for i in 0:(s.set.segments-2)]...)          # connection between tether points
    conn = vcat(conn, [(6+s.set.segments-1, 1)])                                # connection final tether point to bridle point
    return conn  
end   
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Interpolating polars using Dierckx
# ------------------------------------
alpha_cl = [-180.0, -160.0,  -90.0,  -20.0,  -10.0,  -5.0,  0.0, 20.0,  40.0,  90.0, 160.0, 180.0]
cl_list  = [   0.0,    0.5,    0.0,   0.08,  0.125,  0.15,  0.2,  1.0,   1.0,   0.0,  -0.5,   0.0]
alpha_cd = [-180.0, -170.0, -140.0,  -90.0,  -20.0,   0.0, 20.0, 90.0, 140.0, 170.0, 180.0]
cd_list  = [   0.5,    0.5,    0.5,    1.0,    0.2,   0.1,  0.2,  1.0,   0.5,   0.5,   0.5] 
function cl_interp(alpha)
    cl_spline = Spline1D(alpha_cl, cl_list)
    return cl_spline(alpha)
end
function cd_interp(alpha)
    cd_spline = Spline1D(alpha_cd, cd_list)
    return cd_spline(alpha)
end
@register_symbolic cl_interp(alpha)
@register_symbolic cd_interp(alpha)
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Initialize the simulation
# -----------------------------------------
function init_sim!(s::KPS5)
    pos, vel = calc_initial_state(s)
    dt = 1/s.set.sample_freq
    simple_sys,  pos, vel, e_x, e_y, e_z, v_app_point, alpha1p  = model(s, pos, vel)
    s.sys = simple_sys
    tspan = (0.0, s.set.sim_time)
    s.prob = ODEProblem(simple_sys, nothing, tspan)
    #differential_vars = ones(Bool, length(y0))
    #prob = DAEProblem{true}(residual!, yd0, y0, tspan, s; differential_vars)
    #solver  = DFBDF(autodiff=AutoFiniteDiff(), max_order=Val{s.set.max_order}()) 
    #s.integrator = OrdinaryDiffEqCore.init(s.prob, Rodas5(autodiff=false); s.set.dt, abstol=s.set.tol, save_on=false)
    #s.integrator = OrdinaryDiffEqCore.init(prob, solver; abstol=abstol, reltol=s.set.rel_tol, save_everystep=false, initializealg=OrdinaryDiffEqCore.NoInit())
    s.integrator = OrdinaryDiffEqCore.init(s.prob, FBDF(autodiff=false); dt, abstol=s.set.abs_tol, reltol = s.set.rel_tol, save_on=false)
end
# ------------------------------
# Calculate Initial State
# ------------------------------
function calc_initial_state(s::KPS5)  
    p1location = [s.set.l_tether*cos(deg2rad(s.set.elevation)) 0 s.set.l_tether*sin(deg2rad(s.set.elevation))]
    kitepos0rot = get_kite_points(s)
    POS0 = kitepos0rot .+ p1location'
    POS0 = hcat(POS0, zeros(3, 1))
    if s.set.segments > 1
        extra_nodes = [POS0[:,6] + (POS0[:,1] - POS0[:,6]) * i / s.set.segments for i in 1:(s.set.segments-1)]
        POS0 = hcat(POS0, extra_nodes...)
    end     
    VEL0 = zeros(3, points(s))
    return POS0, VEL0
end
# -----------------------------------------------------
# initializing kite points
# -----------------------------------------------------
function get_kite_points(s::KPS5)
    # Original kite points in local reference frame
    kitepos0 =                         # KITE points  
    # P1 Bridle        P2                                    P3                                  P4                              P5
    [0.000         s.set.cord_length/2               -s.set.cord_length/2                       0                                0;
    0.000               0                                    0                                -s.set.width/2               s.set.width/2;
    0.000     s.set.height_k+s.set.h_bridle    s.set.height_k+s.set.h_bridle                 s.set.h_bridle               s.set.h_bridle]

    beta = deg2rad(s.set.elevation)
    Y_r = [sin(beta) 0 cos(beta);
                 0    1       0;
          -cos(beta) 0 sin(beta)]    
    # Apply rotation to all points
    kitepos0rot =  Y_r * kitepos0 
    
    return kitepos0rot
end
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define a function to create reference frame update equations
# ---------------------------------------------------------------
function create_reference_frame_equations(pos, e_x, e_y, e_z)
    # Calculate vectors for the reference frame
    X = pos[:, 2] - pos[:, 3]  # Vector from P3 to P2
    Y = pos[:, 5] - pos[:, 4]  # Vector from P4 to P5
    Z = cross(X, Y)            # Cross product for Z axis
    # Normalize these vectors to get unit vectors
    norm_X = sqrt(sum(X .^ 2))
    norm_Y = sqrt(sum(Y .^ 2))
    norm_Z = sqrt(sum(Z .^ 2))
    # Create equations to update the reference frame
    ref_frame_eqs = [
        e_x[1] ~ X[1] / norm_X,
        e_x[2] ~ X[2] / norm_X,
        e_x[3] ~ X[3] / norm_X,
        e_y[1] ~ Y[1] / norm_Y,
        e_y[2] ~ Y[2] / norm_Y,
        e_y[3] ~ Y[3] / norm_Y,
        e_z[1] ~ Z[1] / norm_Z,
        e_z[2] ~ Z[2] / norm_Z,
        e_z[3] ~ Z[3] / norm_Z
    ]  
    return ref_frame_eqs
end
# --------------------------------------------------------------------------
# computing angle of attack
# --------------------------------------------------------------------------
function compute_alpha1p(v_a, e_z, e_x)
    # Calculate the angle Alpha1p
    v_a_z = v_a ⋅ e_z  # Use the symbolic dot product
    v_a_x = v_a ⋅ e_x  # Use the symbolic dot product
    
    # Use atan for the symbolic computation
    alpha1p =rad2deg.(atan(v_a_z, v_a_x))
    return alpha1p
end
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate Rest Lengths
# ---------------------------
function calc_rest_lengths(s::KPS5)
    conn = getconnections(s)
    POS0, VEL0  = calc_initial_state(s)  
    lengths = [norm(POS0[:,conn[i][2]] - POS0[:,conn[i][1]]) for i in 1:9]
    l10 = norm(POS0[:,1] - POS0[:,6])
    lengths = vcat(lengths, [(l10 + s.set.v_reel_out*t)/s.set.segments for _ in 1:s.set.segments]...)
    return lengths, l10
end
function get_wind_vector(s::KPS5)
    v_wind_magnitude = s.set.v_wind
    wind_angle = deg2rad(s.set.upwind_dir)           # angle measuring form North (pos x), CW(+)
    wind_vector = v_wind_magnitude*[-cos(wind_angle),-sin(wind_angle), 0.0]
    return wind_vector
end
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Define the Model
# -----------------------------------------------
function model(s::KPS5, pos, vel)
    POS0, VEL0 = pos, vel
    rest_lengths, l_tether = calc_rest_lengths(s)
    # @parameters K1=s.set.springconstant_tether K2=s.set.springconstant_bridle K3=s.set.springconstant_kite C1=s.set.damping_tether C2=s.set.rel_damping_bridle*s.set.damping_tether C3=s.set.rel_damping_kite*s.set.damping_tether
    # same as Uwe here, only rel c and k , wrt tether, can/should be improved, to be for kite/bridle separately
    @parameters K1=s.set.c_spring K2=s.set.c_spring K3=s.set.c_spring  C1=s.set.damping C2=s.set.damping*s.set.rel_damping C3=s.set.damping*s.set.rel_damping
    @parameters m_kite=s.set.mass kcu_mass=s.set.kcu_mass rho_tether=s.set.rho_tether rel_compr_k=s.set.rel_compr_stiffness 
    @parameters rho=s.set.rho_0 g_earth=-9.81 cd_tether=s.set.cd_tether d_tether=s.set.d_tether S=s.set.area
    @parameters kcu_cd=s.set.cd_kcu kcu_diameter=s.set.kcu_diameter
    @variables pos(t)[1:3, 1:points(s)] = POS0
    @variables vel(t)[1:3, 1:points(s)] = VEL0
    @variables acc(t)[1:3, 1:points(s)]
    @variables v_app_point(t)[1:3, 1:points(s)]
    @variables segment(t)[1:3, 1:total_segments(s)]
    @variables unit_vector(t)[1:3, 1:total_segments(s)]
    @variables norm1(t)[1:total_segments(s)]
    @variables rel_vel(t)[1:3, 1:total_segments(s)]
    @variables spring_vel(t)[1:total_segments(s)]
    @variables k_spring(t)[1:total_segments(s)]
    @variables spring_force(t)[1:3, 1:total_segments(s)]
    @variables v_apparent(t)[1:3, 1:total_segments(s)]
    @variables v_app_perp(t)[1:3, 1:total_segments(s)]
    @variables norm_v_app(t)[1:total_segments(s)]
    @variables half_drag_force(t)[1:3, 1:total_segments(s)]
    @variables drag_force(t)[1:3, 1:total_segments(s)]
    @variables total_force(t)[1:3, 1:points(s)]
    # local kite reference frame
    @variables e_x(t)[1:3]
    @variables e_y(t)[1:3]
    @variables e_z(t)[1:3] 
    @variables alpha1p(t)[1:1]  

    eqs1 = vcat(D.(pos) .~ vel,
                D.(vel) .~ acc)
    eqs2 = vcat(eqs1...)
    eqs2 = vcat(eqs2, acc[:,6] .~ [0.0, 0.0, 0.0])      # origin is six, make 6 not being hardcoded
    # -----------------------------
    # defining the connections and their respective rest lengths, unit spring constants, damping and masses
    # -----------------------------
                            # connections   adding segment connections, from origin to bridle 
    conn = getconnections(s)
     # final connection last tether point to bridle point
                            # unit spring constants (K1 tether, K2 bridle, K3 kite)
    k_segments = [K2, K3, K2, K2, K3, K3, K2, K3, K3]
    k_segments = vcat(k_segments, [K1 for _ in 1:s.set.segments]...)
                            # unit damping constants (C1 tether, C2 bridle, C3 kite)
    c_segments = [C2, C3, C2, C2, C3, C3, C2, C3, C3]
    c_segments = vcat(c_segments, [C1 for _ in 1:s.set.segments]...)
                            # masses
    mass_bridlelines = ((s.set.d_line/2000)^2)*pi*rho_tether*s.set.l_bridle #total mass entire bridle 
    mass_halfbridleline = mass_bridlelines/8 # half the connection of bridle line to kite (to assign to each kitepoint) so the other 4 halves get assigned to bridlepoint 
    mass_tether = ((d_tether/2000)^2)*pi*rho_tether*l_tether
    mass_tetherpoints = mass_tether/(s.set.segments+1)
    mass_bridlepoint = 4*mass_halfbridleline + kcu_mass + mass_tetherpoints # 4 bridle connections, kcu and tether
    m_kitepoints = (m_kite/4) + mass_halfbridleline 
    PointMasses = [mass_bridlepoint, m_kitepoints, m_kitepoints, m_kitepoints, m_kitepoints]
    PointMasses = vcat(PointMasses, [mass_tetherpoints for _ in 1:s.set.segments]...)
    # getting wind vector
    wind_vector = get_wind_vector(s)
    # -----------------------------
    # Equations for Each Segment (Spring Forces, Drag, etc.)
    # -----------------------------
    for i in 1:total_segments(s)     
        eqs = [
            segment[:, i]      ~ pos[:, conn[i][2]] - pos[:, conn[i][1]],
            norm1[i]           ~ norm(segment[:, i]),
            unit_vector[:, i]  ~ -segment[:, i] / norm1[i],
            rel_vel[:, i]      ~ vel[:, conn[i][2]] - vel[:, conn[i][1]],
            spring_vel[i]      ~ -unit_vector[:, i] ⋅ rel_vel[:, i],
            k_spring[i]        ~ (k_segments[i]/rest_lengths[i]) * (rel_compr_k + (1-rel_compr_k)*(norm1[i] > rest_lengths[i])),
            spring_force[:, i] ~ (k_spring[i]*(norm1[i] - rest_lengths[i]) + c_segments[i] * spring_vel[i]) * unit_vector[:, i],
            v_apparent[:, i]   ~ wind_vector - (vel[:, conn[i][1]] + vel[:, conn[i][2]]) / 2,  # change wind for winddirection! integrate
            v_app_perp[:, i]   ~ v_apparent[:, i] - (v_apparent[:, i] ⋅ unit_vector[:, i]) .* unit_vector[:, i],
            norm_v_app[i]      ~ norm(v_app_perp[:, i])
        ]
        if i > 9 # tether segments
            push!(eqs, half_drag_force[:, i] ~ 0.25 * rho * cd_tether * norm_v_app[i] * (rest_lengths[i]*d_tether/1000) * v_app_perp[:, i])
        elseif i in [1, 3, 4, 7] # bridle lines, try to find Cd_bridlelines later
            push!(eqs, half_drag_force[:, i] ~ 0.25 * rho * cd_tether * norm_v_app[i] * (rest_lengths[i]*(s.set.d_line/1000)) * v_app_perp[:, i])
        else    # kite
            push!(eqs, half_drag_force[:, i] ~ zeros(3))
        end
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
    end
    # -----------------------------
    # Reference Frame and Aerodynamic Coefficients
    # -----------------------------
    ref_frame_eqs = create_reference_frame_equations(pos, e_x, e_y, e_z)
    eqs2 = vcat(eqs2, ref_frame_eqs)
    # only 1 AOA
    v_a_kite = wind_vector - (vel[:, 2] + vel[:, 3])/2     # computing AOA only for center chord   # Appas.set.t wind velocity
    alpha1p = compute_alpha1p(v_a_kite, e_z, e_x)   # Calculate Alpha1p at this time step
    eqs2 = vcat(eqs2, alpha1p[1] ~ alpha1p)  # Add the equation for Alpha1p for each of 4 kite points (first bering bridle so i-1)   
    # getting Cl and Cd
    Cl = cl_interp(alpha1p)            
    Cd = cd_interp(alpha1p)

    # -----------------------------
    # Force Balance at Each Point
    # -----------------------------
    for i in 1:points(s)  
        eqs = []  
        force = sum([spring_force[:, j] for j in 1:total_segments(s) if conn[j][2] == i]; init=zeros(3)) -
                sum([spring_force[:, j] for j in 1:total_segments(s) if conn[j][1] == i]; init=zeros(3)) +
                sum([half_drag_force[:, j] for j in 1:total_segments(s) if conn[j][1] == i]; init=zeros(3)) +
                sum([half_drag_force[:, j] for j in 1:total_segments(s) if conn[j][2] == i]; init=zeros(3))
        v_app_point[:, i] ~ wind_vector - vel[:, i]
        if i == 1                  # KCU drag at bridle point
            area_kcu = pi * ((kcu_diameter / 2) ^ 2)
            Dx_kcu = 0.5*rho*kcu_cd *area_kcu*(v_app_point[1, i]*v_app_point[1, i])
            Dy_kcu = 0.5*rho*kcu_cd *area_kcu*(v_app_point[2, i]*v_app_point[2, i])
            Dz_kcu = 0.5*rho*kcu_cd *area_kcu*(v_app_point[3, i]*v_app_point[3, i])
            D = [Dx_kcu, Dy_kcu, Dz_kcu]
            push!(eqs, total_force[:, i] ~ force + D)
        elseif i in 2:5           # the kite points that get Aero Forces

            v_app_mag_squared = v_app_point[1, i]^2 + v_app_point[2, i]^2 + v_app_point[3, i]^2
            # # Lift calculation
            L_perpoint = (1/4) * 0.5 * rho * Cl * S * (v_app_mag_squared)
            # Cross product and normalization
            cross_vapp_X_e_y = cross(v_app_point[:, i], e_y)
            normcross_vapp_X_e_y = norm(cross_vapp_X_e_y)
            L_direction = cross_vapp_X_e_y / normcross_vapp_X_e_y
             # Final lift force vector
            L = L_perpoint * L_direction

            # Drag calculation
            D_perpoint = (1/4) * 0.5 * rho * Cd * S * v_app_mag_squared
            # # Create drag direction components
            D_direction = [v_app_point[1, i] / norm(v_app_point[:, i]), v_app_point[2, i] / norm(v_app_point[:, i]), v_app_point[3, i] / norm(v_app_point[:, i])]
            # # Final drag force vector components
            D = D_perpoint * D_direction
            
            # Total aerodynamic force
            Fa = [L[1]+ D[1], L[2]+ D[2], L[3]+ D[3]]
            push!(eqs, total_force[1, i] ~ force[1] + Fa[1])
            push!(eqs, total_force[2, i] ~ force[2] + Fa[2])
            push!(eqs, total_force[3, i] ~ force[3] + Fa[3])
        elseif i != 6                      
            push!(eqs, total_force[:, i] ~ force)
        end
        push!(eqs, acc[:, i] ~ [0.0, 0.0, g_earth] + total_force[:, i] / PointMasses[i])
        eqs2 = vcat(eqs2, reduce(vcat, eqs))
        eqs2 = vcat(eqs2, v_app_point[:, i] ~ wind_vector - vel[:, i])
    end
    @named sys = ODESystem(reduce(vcat, Symbolics.scalarize.(eqs2)), t)
    simple_sys = structural_simplify(sys) 
    simple_sys, pos, vel, e_x, e_y, e_z, v_app_point, alpha1p 
end
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# next step function
# -----------------------------
function next_step!(s::KPS5; dt=(1/s.set.sample_freq))
    s.t_0 = s.integrator.t
    steptime = @elapsed OrdinaryDiffEqCore.step!(s.integrator, dt, true)
    s.iter += 1
    s.integrator.t, steptime
end
function simulate(s, logger)
    dt = 1/s.set.sample_freq
    tol = s.set.abs_tol
    tspan = (0.0, dt)
    time_range = 0:dt:s.set.sim_time-dt
    steps = length(time_range)
    iter = 0
    for i in 1:steps
        next_step!(s; dt=dt)
        u = s.get_state(s.integrator)
        x = u[1][1, :]
        y = u[1][2, :]
        z = u[1][3, :]
        iter += s.iter
        sys_state = SysState{points(s)}()
        sys_state.X .= x
        sys_state.Y .= y
        sys_state.Z .= z
        println("iter: $iter", " steps: $steps")
        log!(logger, sys_state)
        println(x[end], ", ", sys_state.X[end])
    end
    println("iter: $iter", " steps: $steps")
    return nothing
end

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Generate Getters
# -----------------------------
function generate_getters!(s::KPS5)
    sys = s.sys
    c = collect
    get_state = ModelingToolkit.getu(sys, 
        [c(sys.pos)]
    )
    s.get_state = (integ) -> get_state(integ)
    return nothing
end


# =================== getter functions ====================================================

"""
    calc_height(s::KPS5)

Determine the height of the topmost kite particle above ground.
"""
function calc_height(s::KPS5)
    pos_kite(s)[3]
end

"""
    pos_kite(s::KPS5)

Return the position of the kite (top particle).
"""
function pos_kite(s::KPS5)
    s.pos[end-2]
end

"""
    kite_ref_frame(s::KPS5; one_point=false)

Returns a tuple of the x, y, and z vectors of the kite reference frame.
"""
function kite_ref_frame(s::KPS5; one_point=false)
    if one_point
        c = s.z
        y = normalize(s.v_apparent × c)
        x = normalize(y × c)
        return x, y, c
    else
        return s.x, s.y, s.z
    end
end

"""
    winch_force(s::KPS5)

Return the absolute value of the force at the winch as calculated during the last timestep. 
"""
function winch_force(s::KPS5) norm(s.last_force) end

"""
    cl_cd(s::KPS5)

Calculate the lift and drag coefficients of the kite, based on the current angles of attack.
"""
function cl_cd(s::KPS5)
    rel_side_area = s.set.rel_side_area/100.0  # defined in percent
    K = 1 - rel_side_area                      # correction factor for the drag
    if s.set.version == 3
        drag_corr = 1.0
    else
        drag_corr = DRAG_CORR
    end
    CL2, CD2 = s.calc_cl(s.alpha_2), drag_corr * s.calc_cd(s.alpha_2)
    CL3, CD3 = s.calc_cl(s.alpha_3), drag_corr * s.calc_cd(s.alpha_3)
    CL4, CD4 = s.calc_cl(s.alpha_4), drag_corr * s.calc_cd(s.alpha_4)
    if s.set.version == 3
        return CL2, CD2
    else
        return CL2, K*(CD2+CD3+CD4)
    end
end

# ==================== end of getter functions ================================================

function spring_forces(s::KPS5; prn=true)
    forces = zeros(SimFloat, s.set.segments+KITE_SPRINGS)
    for i in 1:s.set.segments
        forces[i] =  s.springs[i].c_spring * (norm(s.pos[i+1] - s.pos[i]) - s.segment_length) * s.stiffness_factor
        if forces[i] > s.set.max_force && prn
            println("Tether raptures for segment $i !")
        end
    end
    for i in 1:KITE_SPRINGS
        p1 = s.springs[i+s.set.segments].p1  # First point nr.
        p2 = s.springs[i+s.set.segments].p2  # Second point nr.
        pos1, pos2 = s.pos[p1], s.pos[p2]
        spring = s.springs[i+s.set.segments]
        l_0 = spring.length # Unstressed length
        k = spring.c_spring * s.stiffness_factor       # Spring constant 
        segment = pos1 - pos2
        norm1 = norm(segment)
        k1 = 0.25 * k # compression stiffness kite segments
        if (norm1 - l_0) > 0.0
            spring_force = k *  (norm1 - l_0) 
        else 
            spring_force = k1 *  (norm1 - l_0)
        end
        forces[i+s.set.segments] = spring_force
        if norm(s.spring_force) > 4000.0
            println("Bridle brakes for spring $i !")
        end
    end
    forces
end


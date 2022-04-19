```@meta
CurrentModule = KiteModels
```
## Configuration
To configure the parameters of the kite models, edit the file data/settings.yaml , or create
a copy under a different name and change the name of the active configuration in the file data/system.yaml .

## Parameters
The following parameters are used by this package:
```yaml
system:
    sample_freq: 20        # sample frequency in Hz

initial:
    l_tether: 392.0        # initial tether length       [m]
    elevation: 70.7        # initial elevation angle   [deg]
    v_reel_out: 0.0        # initial reel out speed    [m/s]
    depower:   25.0        # initial depower settings    [%]

solver:
    abs_tol: 0.0006        # absolute tolerance of the DAE solver [m, m/s]
    rel_tol: 0.001         # relative tolerance of the DAE solver [-]
    max_iter:  200         # max number of iterations of the steady-state-solver

steering:
    c0:       0.0          # steering offset   -0.0032           [-]
    c_s:      2.59         # steering coefficient one point model
    c2_cor:   0.93         # correction factor one point model
    k_ds:     1.5          # influence of the depower angle on the steering sensitivity

depower:
    alpha_d_max:    31.0   # max depower angle                            [deg]
    
kite:
    model: "data/kite.obj" # 3D model of the kite
    mass:  6.2             # kite mass incl. sensor unit [kg]
    area: 10.18            # projected kite area         [m²]
    rel_side_area: 30.6    # relative side area           [%]
    height: 2.23           # height of the kite           [m]
    alpha_cl:  [-180.0, -160.0, -90.0, -20.0, -10.0,  -5.0,  0.0, 20.0, 40.0, 90.0, 160.0, 180.0]
    cl_list:   [   0.0,    0.5,   0.0,  0.08, 0.125,  0.15,  0.2,  1.0,  1.0,  0.0,  -0.5,   0.0]
    alpha_cd:  [-180.0, -170.0, -140.0, -90.0, -20.0, 0.0, 20.0, 90.0, 140.0, 170.0, 180.0]
    cd_list:   [   0.5,    0.5,    0.5,   1.0,   0.2, 0.1,  0.2,  1.0,   0.5,   0.5,   0.5]
    
kps4:
    width:         2.23     # width of the kite                      [m]
    alpha_zero:    4.0      # should be 5                      [degrees]
    alpha_ztip:   10.0      #                                  [degrees]
    m_k:           0.2      # relative nose distance; increasing m_k increases C2 of the turn-rate law
    rel_nose_mass: 0.47     # relative nose mass
    rel_top_mass:  0.4      # mass of the top particle relative to the sum of top and side particles

bridle:
    d_line:    2.5         # bridle line diameter                  [mm]
    l_bridle: 33.4         # sum of the lengths of the bridle lines [m]
    h_bridle:  4.9         # height of bridle                       [m]

kcu:
    kcu_mass: 8.4                # mass of the kite control unit   [kg]

tether:
    d_tether:  4           # tether diameter                 [mm]
    cd_tether: 0.958       # drag coefficient of the tether
    damping: 473.0         # unit damping coefficient        [Ns]
    c_spring: 614600.0     # unit spring constant coefficient [N]
    rho_tether:  724.0     # density of Dyneema           [kg/m³]

environment:
    v_wind: 9.51             # wind speed at reference height          [m/s]
    h_ref:  6.0              # reference height for the wind speed     [m]

    rho_0:  1.225            # air density at the ground or zero       [kg/m³]
    alpha:  0.08163          # exponent of the wind profile law
    z0:     0.0002           # surface roughness                       [m]
    profile_law: 3           # 1=EXP, 2=LOG, 3=EXPLOG
```
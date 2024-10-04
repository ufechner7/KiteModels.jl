using KiteModels

kps4::KPS4 = KPS4(KCU(se()))

height = 100.0
v_wind_gnd = kps4.set.v_wind
upwind_dir = -pi - deg2rad(2.0)
wind_dir = -upwind_dir - pi/2
KiteModels.set_v_wind_ground!(kps4, height, v_wind_gnd, wind_dir)
println(rad2deg(KiteModels.upwind_dir(kps4)))

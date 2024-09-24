# test case for calc_aero_forces

function calc_aoa(s::KPS4, pos, vel, alpha_depower=0.0, rel_steering=0.0)
    # pos_B, pos_C, pos_D: position of the kite particles B, C, and D
    # v_B,   v_C,   v_D:   velocity of the kite particles B, C, and D
    pos_B, pos_C, pos_D = pos[s.set.segments+3], pos[s.set.segments+4], pos[s.set.segments+5]
    v_B,   v_C,   v_D   = vel[s.set.segments+3], vel[s.set.segments+4], vel[s.set.segments+5]
    va_2,  va_3,  va_4  = s.v_wind - v_B, s.v_wind - v_C, s.v_wind - v_D
 
    pos_centre = 0.5 * (pos_C + pos_D)
    delta = pos_B - pos_centre
    z = -normalize(delta)
    y = normalize(pos_C - pos_D)
    x = y × z
    s.x .= x; s.y .= y; s.z .= z # save the kite reference frame in the state

    va_xz2 = va_2 - (va_2 ⋅ y) * y
    va_xy3 = va_3 - (va_3 ⋅ z) * z
    va_xy4 = va_4 - (va_4 ⋅ z) * z
    alpha_2 = rad2deg(π - acos2(normalize(va_xz2) ⋅ x) - alpha_depower)     + s.set.alpha_zero
    alpha_3 = rad2deg(π - acos2(normalize(va_xy3) ⋅ x) - rel_steering * s.ks) + s.set.alpha_ztip
    alpha_4 = rad2deg(π - acos2(normalize(va_xy4) ⋅ x) + rel_steering * s.ks) + s.set.alpha_ztip
    alpha_2, alpha_3, alpha_4
end

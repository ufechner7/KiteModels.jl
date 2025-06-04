# SPDX-FileCopyrightText: 2024 Bart van de Lint
# SPDX-License-Identifier: MIT

using LinearAlgebra
using Dierckx
using Plots

# Define constants
rho = 1.0
tip_length = 0.5
middle_length = 1.5 # length of kite in middle
middle_line = 150.0 # length of line
left_line = 150.0 # length of line
right_line = 150.0 # length of line
d_s = 2.0 # distance between closest steering bridles
r = 10.0
w = pi*10.0

# Define vectors
E = [0; 101; 100]
C = [-1; 101.1; 101]
D = [1; 101.1; 101]
v_c = [0; -1; 0]
v_d = [0; -5; 0]
v_e = [0; -1; 0]
v_wind = [0; 10; 0]

alpha_cl = [-180.0, -160.0, -90.0, -20.0, -10.0,  -5.0,  0.0, 20.0, 40.0, 90.0, 160.0, 180.0]
cl_list = [   0.0,    0.5,   0.0,  0.08, 0.125,  0.15,  0.2,  1.0,  1.0,  0.0,  -0.5,   0.0]
alpha_cd = [-10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
cd_list = [   0.2,    0.2,    0.1,   0.1,   0.1, 0.0,  0.1,  0.1,   0.2,   0.2,   0.3]

c_l = Spline1D(deg2rad.(alpha_cl), cl_list)
c_d = Spline1D(deg2rad.(alpha_cd), cd_list)

α_l = π/2 - d_s/(2*r)
α_r = π/2 + d_s/(2*r)

function acos2(arg)
    arg2 = min(max(arg, -one(arg)), one(arg))
    acos(arg2)
 end

function calc_lift(; n=10, left_line=150.0, right_line=150.0)
    # Define functions
    P_c = 0.5.*(C+D)
    e_y = (C - D) ./ norm(C - D)
    e_z = (E - P_c) ./ norm(E - P_c)
    e_x = cross(e_y, e_z)

    F(α) = E + e_y*cos(α)*r - e_z*sin(α)*r
    e_r(α) = (E - F(α))/norm(E-F(α))

    v_cx = dot(v_c, e_x).*e_x
    v_dx = dot(v_d, e_x).*e_x
    v_dy = dot(v_d, e_y).*e_y
    v_dz = dot(v_d, e_z).*e_z
    v_cy = dot(v_c, e_y).*e_y
    v_cz = dot(v_c, e_z).*e_z
    y_lc = norm(C - P_c)
    y_ld = -norm(D - P_c)

    y_l(α) = cos(α) * r
    v_kite(α) = α < π/2 ?
        ((v_cx - v_dx)./(y_lc - y_ld).*(y_l(α) - y_ld) + v_dx) + v_cy + v_cz :
        ((v_cx - v_dx)./(y_lc - y_ld).*(y_l(α) - y_ld) + v_dx) + v_dy + v_dz
    v_a(α) = v_wind - v_kite(α)
    v_a_xr(α) = v_a(α) - (dot(v_a(α), cross(e_r(α), e_x))*cross(e_r(α), e_x))

    length(α) = α < π/2 ?
        (tip_length + (middle_length-tip_length)*α*r/(0.5*w)) :
        (tip_length + (middle_length-tip_length)*(π-α)*r/(0.5*w))

    function d(α)
        if α < α_l
            return middle_line - left_line
        elseif α > α_r
            return middle_line - right_line
        else
            return (-right_line + left_line) / (α_r - α_l) * (α - α_l) + (middle_line - left_line)
        end
    end
    aoa(α) = v_a_xr(α) != [0.0, 0.0, 0.0] ?
        π - acos2(normalize(v_a_xr(α)) ⋅ e_x) + atan(d(α)/length(α)) :
        atan(d(α)/length(α))

    dL_dα(α) = 0.5*rho*(norm(v_a_xr(α)))^2*r*length(α)*c_l(aoa(α)) .* normalize(v_a_xr(α) × (e_r(α) × e_x))
    dD_dα(α) = 0.5*rho*norm(v_a_xr(α))*r*length(α)*c_d(aoa(α)) .* v_a_xr(α) # the sideways drag cannot be calculated with the C_d formula

    # Calculate the integral
    α_0 = pi/2 - w/2/r
    α_middle = pi/2
    dα = (α_middle - α_0) / n
    L_C = sum(dL_dα(α_0 + dα/2 + i*dα) * dα for i in 1:n)
    L_D = sum(dL_dα(pi - (α_0 + dα/2 + i*dα)) * dα for i in 1:n)
    D_C = sum(dD_dα(α_0 + dα/2 + i*dα) * dα for i in 1:n)
    D_D = sum(dD_dα(pi - (α_0 + dα/2 + i*dα)) * dα for i in 1:n)
    println(L_C, L_D, D_C, D_D);
    return plot(1:n*2, norm.(dL_dα(α_0 + dα/2 + i*dα) * dα for i in 1:n*2), title="Lift per segment")
end

# println(calc_lift())
# println(calc_lift(n=10, left_line=149, right_line=149.9))
savefig(calc_lift(n=100, left_line=149.9, right_line=149), "n=100")
savefig(calc_lift(n=5, left_line=149.9, right_line=149), "n=5")
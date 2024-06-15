using LinearAlgebra
using QuadGK
using Dierckx

# Define constants
rho = 1.0
tip_length = 0.5
middle_length = 1.5 # length of kite in middle
middle_line = 150.0 # length of line
left_line = 150.0 # length of line
right_line = 150.0 # length of line
d_s = 1.0 # distance between closest steering bridles
r = 10.0
w = pi*10.0

# Define vectors
E = [0; 101; 100]
C = [-1; 101; 101]
D = [1; 101; 101]
v_c = [0; -1; 0]
v_d = [0; -1; 0]
v_e = [0; -1; 0]
v_wind = [0; 1; 0]



const delta_cl = [-10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0] # steering_line_length - middle_line_length in cm
const cl_list = [0.15, 1.0, 1.0, 0.08, 0.125, 0.15, 0.0, -0.1, -0.5, -0.5, -1.0]
const c_l = Spline1D(delta_cl./100, cl_list)

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
    v_ey = dot(v_e, e_y).*e_y
    v_ez = dot(v_e, e_z).*e_z
    y_lc = norm(C - P_c)
    y_ld = -norm(D - P_c)

    y_l(α) = cos(α) * r
    v_kite(α) = (v_cx - v_dx)./(y_lc - y_ld).*(y_l(α) - y_ld) + v_dx + v_ey + v_ez
    v_a(α) = v_wind - v_kite(α)
    v_a_xr(α) = v_a(α) - (dot(v_a(α), cross(e_r(α), e_z))*cross(e_r(α), e_z))
    length(α) = (tip_length + (middle_length-tip_length)*α*r/(0.5*w))

    α_l = 0.5*pi - d_s/(2*r)
    α_r = 0.5*pi + d_s/(2*r)
    function d(α)
        if α < α_l
            println(middle_line - left_line)
            return middle_line - left_line
        elseif α > α_r
            return middle_line - right_line
        else
            return (left_line - right_line) / (α_r - α_l) * (α - α_l) + (middle_line - right_line)
        end
    end

    dL_dα(α) = (0.5*rho*(norm(v_a_xr(α)))^2*r*length(α)*c_l(d(α))*e_r(α))


    # Calculate the integral
    # L, err = quadgk(dL_dalpha, 0, pi/2, rtol=1e-1, atol=1e-1)
    a = 0
    b = pi/2
    h = (b - a) / n

    L_C = h/2 * (dL_dα(a) + 2*sum(dL_dα(a + i*h) for i in 1:n-1) + dL_dα(b))
end
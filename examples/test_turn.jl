# test the turn function
UPWIND_DIR = -pi/2 +deg2rad(10)

# function turn(res, upwind_dir)
#     turnangle = upwind_dir + pi/2
#     res2 = zeros(SimFloat, length(res))
#     for i in 1:div(length(res), 3)
#         x = res[3*(i-1)+1]
#         y = res[3*(i-1)+2]
#         res2[3*(i-1)+1] = cos(turnangle) * x + sin(turnangle) * y
#         res2[3*(i-1)+2] = cos(turnangle) * y - sin(turnangle) * x
#         res2[3*(i-1)+3] = res[3*(i-1)+3]
#     end
#     res2
# end

# y, yd
res=[1.0, 0, 1, 0, 0, 0]
res2 = KiteModels.turn(res, UPWIND_DIR)
println("res2: ", res2)
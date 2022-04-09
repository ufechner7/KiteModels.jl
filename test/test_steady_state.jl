using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils, Plots
using KiteModels, KitePodModels

if ! @isdefined kcu
    const kcu = KCU()
end
if ! @isdefined kps4
    const kps4 = KPS4(kcu)
end

const dt = 0.05

function plot2d(x, z; zoom=1)
    if zoom ==1
        x_max=maximum(x)
        z_max=maximum(z)
        plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false, xlims = (x_max-15.0, x_max+5), ylims = (z_max-15.0, z_max+5))
        plot!([x[7],x[10]],[z[7],z[10]], legend=false)
        plot!([x[8],x[11]],[z[8],z[11]], legend=false)
    else
        plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false)
    end
    plot!(x, z, seriestype = :scatter)
end

function simulate(integrator, steps)
    for i in 1:steps
        # println("lift, drag    [N]  : $(KiteModels.lift_drag(kps4))")
        KiteModels.next_step(kps4, integrator, dt)
        if kps4.stiffness_factor < 1.0
            kps4.stiffness_factor+=0.01
        end
        # local x = Float64[] 
        # local z = Float64[]
        # for i in 1:length(kps4.pos)
        #     push!(x, kps4.pos[i][1])
        #     push!(z, kps4.pos[i][3])
        # end
        # local p = plot2d(x, z; zoom=0)
        # display(p)
        # sleep(0.1)
    end
end

integrator=KiteModels.init_sim(kps4, 1.0)
oldpos = deepcopy(kps4.pos)
kps4.stiffness_factor=0.02
simulate(integrator, 3)
@time simulate(integrator, 100)
x = Float64[] 
z = Float64[]
for i in 1:length(kps4.pos)
     push!(x, kps4.pos[i][1])
     push!(z, kps4.pos[i][3])
end
# p = plot2d(x, z; zoom=0)
# display(p)
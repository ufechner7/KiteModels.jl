using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
    pkg"add ControlPlots#main"
end
using ControlPlots
using KiteModels, KitePodModels, KiteUtils

set = deepcopy(se())
kcu::KCU = KCU(set)
kps3::KPS3 = KPS3(kcu)

function plot2d(pos, reltime=0.0; zoom=true, front=false, segments=6)
    x = Float64[] 
    z = Float64[]
    for i in eachindex(pos)
        if front
            push!(x, pos[i][2])
        else
            push!(x, pos[i][1])
        end
        push!(z, pos[i][3])
    end
    x_max = maximum(x)
    z_max = maximum(z)

    xlabel = "x [m]"
    ann = (x_max-10.0, z_max-3.0, "t=$(round(reltime,digits=1)) s")
    ControlPlots.plot(x, z; xlabel, ylabel="z [m]", xlims = (x_max-15.0, x_max+5), ylims = (z_max-15.0, z_max+5), ann, scatter=true)

    # if length(pos) > segments+1
    #     s=segments
    #     plot!([x[s+1],x[s+4]],[z[s+1],z[s+4]], legend=false) # S6
    #     plot!([x[s+2],x[s+5]],[z[s+2],z[s+5]], legend=false) # S8
    #     plot!([x[s+3],x[s+5]],[z[s+3],z[s+5]], legend=false) # S7
    #     plot!([x[s+2],x[s+4]],[z[s+2],z[s+4]], legend=false) # S2
    #     plot!([x[s+1],x[s+5]] ,[z[s+1],z[s+5]],legend=false) # S5
    # end
end

reltime = 0.0
integrator = KiteModels.init_sim!(kps3, stiffness_factor=0.04)
p = plot2d(kps3.pos, reltime; zoom=false, front=false, segments=set.segments)
display(p)   
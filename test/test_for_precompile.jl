let
    using KiteModels, KitePodModels, KiteUtils
    kcu = KCU(se())
    kps4 = KPS4(kcu)
    dt = 0.05
    STATISTIC = false
    FRONT_VIEW = false
    ZOOM = false
    PLOT = true

    if PLOT
        using Plots
        function plot2d(pos, reltime=0.0; zoom=true, front=false, segments=se())
            x = Float64[] 
            z = Float64[]
            for i in 1:length(pos)
                if front
                    push!(x, pos[i][2])
                else
                    push!(x, pos[i][1])
                end
                push!(z, pos[i][3])
            end
            x_max = maximum(x)
            z_max = maximum(z)
            if zoom
                xlabel = "x [m]"
                if front xlabel = "y [m]" end
                plot(x,z, xlabel=xlabel, ylabel="z [m]", legend=false, xlims = (x_max-15.0, x_max+5), ylims = (z_max-15.0, z_max+5))
                annotate!(x_max-10.0, z_max-3.0, "t=$(round(reltime,digits=1)) s")
            else
                plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false)
                annotate!(x_max-10.0, z_max-3.0, "t=$(round(reltime,digits=1)) s")
            end
            if length(pos) > segments+1
                s=segments
                plot!([x[s+1],x[s+4]],[z[s+1],z[s+4]], legend=false) # S6
                plot!([x[s+2],x[s+5]],[z[s+2],z[s+5]], legend=false) # S8
                plot!([x[s+3],x[s+5]],[z[s+3],z[s+5]], legend=false) # S7
                plot!([x[s+2],x[s+4]],[z[s+2],z[s+4]], legend=false) # S2
                plot!([x[s+1],x[s+5]] ,[z[s+1],z[s+5]],legend=false) # S5
            end
            plot!(x, z, seriestype = :scatter) 
        end
    end

    function simulate(integrator, steps, plot=false)
        start = integrator.p.iter
        for i in 1:steps
            lift, drag = KiteModels.lift_drag(kps4)
            KiteModels.next_step!(kps4, integrator)
            if kps4.stiffness_factor < 1.0
                kps4.stiffness_factor+=0.01
            end
            if plot
                reltime = i*dt
                p = plot2d(kps4.pos, reltime; zoom=ZOOM, front=FRONT_VIEW)
                display(p)
            end
        end
        (integrator.p.iter - start) / steps
    end
    integrator = KiteModels.init_sim!(kps4, prn=STATISTIC)
    kps4.stiffness_factor = 0.04
    simulate(integrator, 100, true)
end

@info "Precompile script has completed execution."
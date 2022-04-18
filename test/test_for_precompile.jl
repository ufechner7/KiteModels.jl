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
        include("plot2d.jl")
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
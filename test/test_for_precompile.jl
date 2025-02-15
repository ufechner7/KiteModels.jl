let
    using KiteModels, KitePodModels, KiteUtils
    kcu = KCU(se())
    kps4 = KPS4(kcu)
    dt = 0.05
    STATISTIC = false
    FRONT_VIEW = false
    ZOOM = false
    PLOT = false
    
    function simulate(integrator, steps, plot=false)
        start = integrator.p.iter
        for i in 1:steps  
            KiteModels.next_step!(kps4, integrator; set_speed=0, dt=dt)      
        end
        (integrator.p.iter - start) / steps
    end
    integrator = KiteModels.init_sim!(kps4, prn=STATISTIC)
    kps4.stiffness_factor = 0.04
    simulate(integrator, 100, true)
end

@info "Precompile script has completed execution."
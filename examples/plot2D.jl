# is included from simulate.jl
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
    if front xlabel = "y [m]" end
    if zoom
        ann = (x_max-10.0, z_max-3.0, "t=$(round(reltime,digits=1)) s")
        p=ControlPlots.plot(x, z; xlabel, ylabel="z [m]", xlims = (x_max-15.0, x_max+5), ylims = (z_max-15.0, z_max+5), ann, scatter=true)
    else
        ann = (x_max-10.0, z_max-3.0, "t=$(round(reltime,digits=1)) s")
        p=ControlPlots.plot(x, z; xlabel, ylabel="z [m]", ann, scatter=true)
    end
    display(p)
    if length(pos) > segments+1
        s=segments
        # plot!([x[s+1],x[s+4]],[z[s+1],z[s+4]], legend=false) # S6
        # plot!([x[s+2],x[s+5]],[z[s+2],z[s+5]], legend=false) # S8
        # plot!([x[s+3],x[s+5]],[z[s+3],z[s+5]], legend=false) # S7
        # plot!([x[s+2],x[s+4]],[z[s+2],z[s+4]], legend=false) # S2
        # plot!([x[s+1],x[s+5]] ,[z[s+1],z[s+5]],legend=false) # S5
    end
    sleep(0.001)
end

function plot2d_(pos, reltime=0.0; zoom=true, front=false, segments=6, line, sc, txt)
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
    if front xlabel = "y [m]" end
    if isnothing(line)
        line, = plt.plot(x,z; linewidth="1")
        sc  = plt.scatter(x, z; s=25, color="red") 
        txt = plt.annotate("t=$(round(reltime,digits=1)) s",  
                            xy=(x_max, z_max+8.0), fontsize = 14)
        plt.xlim(0, x_max+20)
        plt.ylim(0, z_max+20)
        plt.xlabel(xlabel, fontsize=14)
        plt.ylabel("z [m]", fontsize=14)
        plt.grid(true)
        plt.grid(which="major", color="#DDDDDD")
    else
        line.set_xdata(x)
        line.set_ydata(z)
        sc.set_offsets(hcat(x,z))
        txt.set_text("t=$(round(reltime,digits=1)) s")
        plt.gcf().canvas.draw()
    end
    sleep(0.01)
    line, sc, txt
end

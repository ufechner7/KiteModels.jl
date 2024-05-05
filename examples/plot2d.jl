function plot2d_(pos, reltime; zoom=true, front=false, segments=6, fig="", lines, sc, txt)
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
    if isnothing(lines)
        if fig != ""
            plt.figure(fig)
        end
        lines=[]
        line, = plt.plot(x,z; linewidth="1")
        push!(lines, line)
        sc  = plt.scatter(x, z; s=25, color="red") 
        if zoom
            txt = plt.annotate("t=$(round(reltime,digits=1)) s",  
                xy=(x_max, z_max+4.7), fontsize = 14)
            plt.xlim(x_max-15.0, x_max+20)
            plt.ylim(z_max-15.0, z_max+8)
        else
            txt = plt.annotate("t=$(round(reltime,digits=1)) s",  
            xy=(x_max, z_max+8.0), fontsize = 14)
            plt.xlim(0, x_max+20)
            plt.ylim(0, z_max+20)
        end
        if length(pos) > segments+1
            s=segments
            line, = plt.plot([x[s+1],x[s+4]],[z[s+1],z[s+4]], linewidth="1"); push!(lines, line) # S6
            line, = plt.plot([x[s+2],x[s+5]],[z[s+2],z[s+5]], linewidth="1"); push!(lines, line) # S8
            line, = plt.plot([x[s+3],x[s+5]],[z[s+3],z[s+5]], linewidth="1"); push!(lines, line) # S7
            line, = plt.plot([x[s+2],x[s+4]],[z[s+2],z[s+4]], linewidth="1"); push!(lines, line) # S2
            line, = plt.plot([x[s+1],x[s+5]],[z[s+1],z[s+5]], linewidth="1"); push!(lines, line) # S5
        end
        plt.xlabel(xlabel, fontsize=14)
        plt.ylabel("z [m]", fontsize=14)
        plt.grid(true)
        plt.grid(which="major", color="#DDDDDD")
    else
        lines[1].set_xdata(x)
        lines[1].set_ydata(z)
        if length(pos) > segments+1
            s=segments
            lines[2].set_xdata([x[s+1],x[s+4]]) # S6
            lines[2].set_ydata([z[s+1],z[s+4]]) # S6
            lines[3].set_xdata([x[s+2],x[s+5]]) # S8
            lines[3].set_ydata([z[s+2],z[s+5]]) # S8
            lines[4].set_xdata([x[s+3],x[s+5]]) # S7
            lines[4].set_ydata([z[s+3],z[s+5]]) # S7
            lines[5].set_xdata([x[s+2],x[s+4]]) # S2
            lines[5].set_ydata([z[s+2],z[s+4]]) # S2
            lines[6].set_xdata([x[s+1],x[s+5]]) # S5
            lines[6].set_ydata([z[s+1],z[s+5]]) # S5
        end
        sc.set_offsets(hcat(x,z))
        txt.set_text("t=$(round(reltime,digits=1)) s")
        plt.gcf().canvas.draw()
    end
    sleep(0.01)
    lines, sc, txt
end

plot2d = let lines = nothing, sc = nothing, txt = nothing  # Note: Must all be on same line as let!
    function(pos, reltime=0.0; kwargs...)
        lines, sc, txt = plot2d_(pos, reltime; lines, sc, txt, kwargs...)
    end
end

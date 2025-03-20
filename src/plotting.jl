
function plot(sys::PointMassSystem, reltime; zoom=false, front = false)
    pos = [sys.points[i].pos_w for i in eachindex(sys.points)]
    seg = [[sys.segments[i].points[1], sys.segments[i].points[2]] for i in eachindex(sys.segments)]
    if zoom && !front
        xlim = (pos[1][1] - 5, pos[1][1]+5)
        ylim = (pos[1][3] - 8, pos[1][3]+2)
    elseif zoom && front
        xlim = (pos[1][2] - 5, pos[1][2]+5)
        ylim = (pos[1][3] - 8, pos[1][3]+2)
    elseif !zoom && !front
        xlim = (0, 60)
        ylim = (0, 60)
    elseif !zoom && front
        xlim = (-30, 30)
        ylim = (0, 60)
    end
    ControlPlots.plot2d(pos, seg, reltime; zoom, front, xlim, ylim)
    # ControlPlots.plot2d(pos, seg, reltime; zoom=false, front=false)
end

function plot(s::KPSQ, reltime; kwargs...)
    pos = s.integrator[s.simple_sys.pos]
    for (i, point) in enumerate(s.point_system.points)
        point.pos_w .= pos[:, i]
    end
    plot(s.point_system, reltime; kwargs...)
end
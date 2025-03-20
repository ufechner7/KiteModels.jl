
function plot(sys::PointMassSystem, reltime)
    pos = [sys.points[i].pos_w for i in eachindex(sys.points)]
    seg = [[sys.segments[i].points[1], sys.segments[i].points[2]] for i in eachindex(sys.segments)]
    ControlPlots.plot2d(pos, seg, reltime; zoom=false, front=true, xlim=(-30, 30), ylim=(0, 60))
    # ControlPlots.plot2d(pos, seg, reltime; zoom=false, front=false)
end

function plot(s::KPSQ, reltime)
    pos = s.integrator[s.simple_sys.pos]
    for (i, point) in enumerate(s.point_system.points)
        point.pos_w .= pos[:, i]
    end
    plot(s.point_system, reltime)
end
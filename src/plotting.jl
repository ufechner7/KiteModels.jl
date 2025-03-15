
function plot(sys::PointMassSystem, reltime)
    pos = [sys.points[i].pos_w for i in eachindex(sys.points)]
    seg = [[sys.segments[i].points[1], sys.segments[i].points[2]] for i in eachindex(sys.segments)]
    ControlPlots.plot2d(pos, seg, reltime; zoom=false, front=true, xlim=(-2, 2), ylim=(0, 120))
    # ControlPlots.plot2d(pos, seg, reltime; zoom=false, front=true)
end

function plot(s::KPSQ, reltime)
    pos = s.integrator[s.simple_sys.pos]
    for (i, point) in enumerate(s.point_system.points)
        point.pos_w .= pos[:, i]
    end
    plot(s.point_system, reltime)
end
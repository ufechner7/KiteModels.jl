
function plot(sys::PointMassSystem, reltime)
    pos = [sys.points[i].pos_w for i in eachindex(sys.points)]
    seg = [[sys.segments[i].points[1], sys.segments[i].points[2]] for i in eachindex(sys.segments)]
    ControlPlots.plot2d(pos, seg, reltime; zoom=false, front=false, xlim=(6, 10), ylim=(45, 60))
    nothing
end
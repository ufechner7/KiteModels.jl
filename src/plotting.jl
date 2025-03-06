
function plot(sys::PointMassSystem, reltime)
    pos = [sys.points[i].pos for i in eachindex(sys.points)]
    seg = [[sys.segments[i].points[1], sys.segments[i].points[2]] for i in eachindex(sys.segments)]
    ControlPlots.plot2d(pos, seg, reltime; zoom=false, front=true, xlim=(-5, 5), ylim=(-20, 5))
    nothing
end
# Copyright (c) 2022, 2024 Uwe Fechner
# SPDX-License-Identifier: MIT

function plot(sys::PointMassSystem, reltime; kite_pos=nothing, zoom=false, front=false)
    pos = [sys.points[i].pos_w for i in eachindex(sys.points)]
    !isnothing(kite_pos) && (pos = [pos..., kite_pos])
    seg = [[sys.segments[i].points[1], sys.segments[i].points[2]] for i in eachindex(sys.segments)]
    if zoom && !front
        xlim = (pos[end][1] - 5, pos[end][1]+5)
        ylim = (pos[end][3] - 8, pos[end][3]+2)
    elseif zoom && front
        xlim = (pos[end][2] - 5, pos[end][2]+5)
        ylim = (pos[end][3] - 8, pos[end][3]+2)
    elseif !zoom && !front
        xlim = (0, 60)
        ylim = (0, 60)
    elseif !zoom && front
        xlim = (-30, 30)
        ylim = (0, 60)
    end
    ControlPlots.plot2d(pos, seg, reltime; zoom, front, xlim, ylim, dz_zoom=0.6)
end

function plot(s::RamAirKite, reltime; kwargs...)
    pos = s.integrator[s.sys.pos]
    kite_pos = s.integrator[s.sys.kite_pos]
    for (i, point) in enumerate(s.point_system.points)
        point.pos_w .= pos[:, i]
    end
    plot(s.point_system, reltime; kite_pos, kwargs...)
end
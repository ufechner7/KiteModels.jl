# SPDX-FileCopyrightText: 2025 Bart van de Lint, Uwe Fechner
#
# SPDX-License-Identifier: MIT

module KiteModelsControlPlotsExt
using ControlPlots, KiteModels

export plot

function ControlPlots.plot(sys::SystemStructure, reltime; l_tether=50.0, wing_pos=nothing, e_z=zeros(3), zoom=false, front=false)
    pos = [sys.points[i].pos_w for i in eachindex(sys.points)]
    !isnothing(wing_pos) && (pos = [pos..., wing_pos...])
    seg = [[sys.segments[i].point_idxs[1], sys.segments[i].point_idxs[2]] for i in eachindex(sys.segments)]
    if zoom && !front
        xlim = (pos[end][1] - 6, pos[end][1]+6)
        ylim = (pos[end][3] - 10, pos[end][3]+2)
    elseif zoom && front
        xlim = (pos[end][2] - 6, pos[end][2]+6)
        ylim = (pos[end][3] - 10, pos[end][3]+2)
    elseif !zoom && !front
        xlim = (-5, l_tether+10)
        ylim = (-5, l_tether+10)
    elseif !zoom && front
        xlim = (-12.5 - 0.5l_tether, 12.5 + 0.5l_tether)
        ylim = (-5, l_tether+10)
    end
    ControlPlots.plot2d(pos, seg, reltime; zoom, front, xlim, ylim, dz_zoom=0.6)
end

function ControlPlots.plot(s::SymbolicAWEModel, reltime; kwargs...)
    wings = s.sys_struct.wings
    pos = s.integrator[s.sys.pos]
    if length(wings) > 0
        wing_pos = [s.integrator[s.sys.wing_pos[i, :]] for i in eachindex(wings)]
        e_z = [s.integrator[s.sys.e_z[i, :]] for i in eachindex(wings)]
    else
        wing_pos = nothing
        e_z = zeros(3)
    end
        
    for (i, point) in enumerate(s.sys_struct.points)
        point.pos_w .= pos[:, i]
    end
    plot(s.sys_struct, reltime; s.set.l_tether, wing_pos, e_z, kwargs...)
end

end

function plot_front_view3(lg)
    display(plotxy(lg.y, lg.z;
    xlabel="pos_y [m]",
    ylabel="height [m]",
    fig="front_view"))
    nothing
end

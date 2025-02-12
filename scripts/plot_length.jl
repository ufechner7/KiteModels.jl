using Serialization

# Plot the interpolations
function plot_interpolations(itp_max, itp_min, gammas, max_xs, min_xs)
    # Generate finer gamma values for smooth interpolation plots
    gamma_fine = range(minimum(gammas), maximum(gammas), length=1000)
    
    # Evaluate interpolations at finer gamma values
    max_x_fine = itp_max.(gamma_fine)
    min_x_fine = itp_min.(gamma_fine)
    
    # Create the plot
    plot(gamma_fine, max_x_fine, label="Maximum x", xlabel="Gamma (radians)", ylabel="x coordinate", title="Interpolated Max and Min x Coordinates", linewidth=2)
    plot!(gamma_fine, min_x_fine, label="Minimum x", linewidth=2)
    
    # Overlay the original data points
    scatter!(gammas, max_xs, label="", markercolor=:blue, markersize=3)
    scatter!(gammas, min_xs, label="", markercolor=:red, markersize=3)
    
    display(plot!())
end

(itp_max, itp_min, gammas, max_xs, min_xs) = deserialize(joinpath(dirname(get_data_path()), se("system_3l.yaml").polar_file))
plot_interpolations(itp_max, itp_min, gammas, max_xs, min_xs)
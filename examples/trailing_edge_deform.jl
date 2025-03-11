"""
Plot the trailing edge deformation
"""

using Pkg, Statistics
using KiteUtils

if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
else
    using ControlPlots
end

se::Settings = KiteUtils.se("system_3l.yaml")

function normalize!(x, y)
    x_min = minimum(x)
    x_max = maximum(x)
    for i in eachindex(x)
        x[i] = (x[i] - x_min) / (x_max - x_min)
        y[i] = (y[i] - x_min) / (x_max - x_min)
    end
end

function turn_flap!(angle, x, y, lower_turn, upper_turn)
    theta = deg2rad(angle)
    x_turn = 0.75
    turn_distance = upper_turn - lower_turn
    smooth_idx = []
    rm_idx = []

    sign = theta > 0 ? 1 : -1
    y_turn = theta > 0 ? upper_turn : lower_turn
    for i in eachindex(x)
        if x_turn - turn_distance < x[i] < x_turn + turn_distance && sign * y[i] > 0
            append!(smooth_idx, i)
        elseif sign * y[i] < 0 && x_turn > x[i] > x_turn - turn_distance
            append!(rm_idx, i)
        end
        if x[i] > x_turn
            x_rel = x[i] - x_turn
            y_rel = y[i] - y_turn
            x[i] = x_turn + x_rel * cos(theta) + y_rel * sin(theta)
            y[i] = y_turn - x_rel * sin(theta) + y_rel * cos(theta)
            if theta > 0 && x[i] < x_turn - turn_distance/2 && y[i] > lower_turn
                append!(rm_idx, i)
            elseif theta < 0 && x[i] < x_turn - turn_distance/2 && y[i] < upper_turn
                append!(rm_idx, i)
            end
        end
    end

    #TODO: lower and upper is slightly off because of smoothing
    lower_i, upper_i = minimum(smooth_idx), maximum(smooth_idx)
    for i in smooth_idx
        window = min(i - lower_i + 1, upper_i - i + 1)
        x[i] = mean(x[i-window:i+window])
    end
    deleteat!(x, rm_idx)
    deleteat!(y, rm_idx)
    nothing
end

function get_lower_upper(x, y)
    lower_flap = 0.0
    upper_flap = 0.0
    min_lower_distance = Inf
    min_upper_distance = Inf
    for (xi, yi) in zip(x, y)
        if yi < 0
            lower_distance = abs(xi - 0.75)
            if lower_distance < min_lower_distance
                min_lower_distance = lower_distance
                lower_flap = yi
            end
        else
            upper_distance = abs(xi - 0.75)
            if upper_distance < min_upper_distance
                min_upper_distance = upper_distance
                upper_flap = yi
            end
        end
    end
    return lower_flap, upper_flap
end

# Select one angle of attack for demonstration
alpha = 0.0

# Create a range of flap angles to visualize
flap_angles = -30:15:30

# Create a plot
rcParams = plt.PyDict(plt.matplotlib."rcParams")
rcParams["text.usetex"] = true

# Read and normalize coordinates
x, y = open(se.foil_file, "r") do f
    x = Float64[]; y = Float64[]
    for line in eachline(f)
        entries = split(chomp(line))
        try
            push!(x, parse(Float64, entries[1]))
            push!(y, parse(Float64, entries[2]))
        catch ArgumentError
        end
    end
    x, y
end
normalize!(x, y)
lower, upper = get_lower_upper(x, y)

# Plot each flap configuration
for angle in flap_angles
    x_mod = copy(x)
    y_mod = copy(y)
    turn_flap!(angle, x_mod, y_mod, lower, upper)
    plt.plot(x_mod, y_mod, label="\$\\alpha\$: $angle", linewidth=1)
end
plt.plot([0.75], [upper], "ro", label="\$y_{\\mathrm{up}}\$", markersize=5)
plt.plot([0.75], [lower], "bo", label="\$y_{\\mathrm{down}}\$", markersize=5)

plt.axis("equal")
plt.legend(loc="upper left")
plt.xlabel("\$x\$")
plt.ylabel("\$y\$")
plt.savefig("test.svg")
plt.show()
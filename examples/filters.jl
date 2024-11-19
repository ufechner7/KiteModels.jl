"""
    ema_filter(measurement, last_measurement, cut_off_freq, dt)

An exponential moving average (EMA), also known as an exponentially weighted moving average (EWMA),
is a type of infinite impulse response filter that applies weighting factors which decrease exponentially.
The weighting for each older datum decreases exponentially, never reaching zero.

# Arguments
- `measurement`: The measurement value
- `last_measurement`: The last measurement value
- `cut_off_freq`: The cut off frequency in Hz
- `dt`: The sampling time in seconds

# Returns
- The filtered value
"""
function ema_filter(measurement, last_measurement, cut_off_freq, dt)
    if cut_off_freq > 0.0
        alpha = dt / (dt + one(measurement) / (2π * cut_off_freq))
        filtered_value = alpha * measurement + (one(measurement) - alpha) * last_measurement
    else
        filtered_value = measurement
    end
    return filtered_value
end
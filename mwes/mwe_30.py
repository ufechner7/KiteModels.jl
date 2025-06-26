import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider

alpha_cl = [-180.0, -160.0, -90.0, -20.0, -10.0, -5.0, 0.0, 5.0, 7.0, 8.2, 12.0, 16.0, 40.0, 90.0, 160.0, 180.0]
cl_list = [0.0, 0.5, 0.0, 0.08, -0.125, 0.05, 0.26, 0.54, 0.68, 0.76, 1.03, 1.13, 0.8, 0.0, -0.5, 0.0]
alpha_cd = [-180.0, -170.0, -140.0, -90.0, -20.0, 0.0, 7.5, 9.5, 12.0, 16.0, 90.0, 140.0, 170.0, 180.0]
cd_list = [0.5, 0.5, 0.5, 1.0, 0.06, 0.06, 0.137, 0.158, 0.19, 0.23, 1.0, 0.5, 0.5, 0.5]

test_alpha = np.linspace(-180, 180, endpoint=True)
test_alpha_2 = np.linspace(-180, 180, len(test_alpha) * 2, endpoint=True)

def update_cl_fit(val):
    return np.poly1d(np.polyfit(test_alpha, np.interp(test_alpha, alpha_cl, cl_list), val))(test_alpha_2)

def update_cd_fit(val):
    return np.poly1d(np.polyfit(test_alpha, np.interp(test_alpha, alpha_cd, cd_list), val))(test_alpha_2)


fig, (ax0, ax1, ax2) = plt.subplots(3, 1)
ax0.plot(alpha_cl, cl_list, 'r-', label='cl interp')
ax0.plot(alpha_cd, cd_list, 'b-', label='cd interp')

line_cl, = ax0.plot(test_alpha_2, update_cl_fit(1), 'r--', label='cl fit')
line_cd, = ax0.plot(test_alpha_2, update_cl_fit(1), 'b--', label='cd fit')

cl_slider = Slider(
    ax=ax1,
    label='n coefficients Cl [-]',
    valmin=1,
    valmax=len(test_alpha),
    valinit=1,
    valstep=1
)
cl_slider.on_changed(lambda n: line_cl.set_ydata(update_cl_fit(n)))

cd_slider = Slider(
    ax=ax2,
    label='n coefficients Cd [-]',
    valmin=1,
    valmax=len(test_alpha),
    valinit=1,
    valstep=1
)
cd_slider.on_changed(lambda n: line_cd.set_ydata(update_cd_fit(n)))

plt.tight_layout()
ax0.legend()

plt.show()

exit(0)

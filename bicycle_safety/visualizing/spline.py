import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
f = CubicSpline([-2, -1, 0, 1, 2], [-2, -1, 0, 10, 2], extrapolate=True)
t = np.linspace(-2, 3)
plt.plot(t, f(t))
plt.show()
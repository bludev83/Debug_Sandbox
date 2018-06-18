"""docstring"""
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 2 * np.pi, 20)
y = np.sin(x)
plt.plot(x, np.sin(x))
plt.show()

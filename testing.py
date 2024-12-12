import numpy as np
import matplotlib.pyplot as plt

def disp(om_t, a, b, om_c):
    return 1 / (om_t - a) ** 2 + 1 / (om_c * (om_t - a * b)) ** 2



a = 1
b = -0.5
om_c = 12

instability_criterion = a * np.abs(b - 1) < (1 + om_c ** (2 / 3)) ** (3 / 2) / om_c

# om_t = np.linspace(-2, 2.5, 1000)

# plt.plot(om_t, disp(om_t, a, b, om_c), label=r'$f(\tilde{\omega})$')
# plt.vlines(a, 0, 100, color='red', linestyle='--', label=r'$\tilde{\omega}=\alpha$')
# plt.vlines(a * b, 0, 100, color='blue', linestyle='--', label=r'$\tilde{\omega}=\alpha\cdot\beta$')
# plt.hlines(1, -2, 2.5, color='black')
# plt.ylim(0, 10)
# plt.xlabel(r'$\tilde{\omega}$')
# plt.ylabel(r'$f(\tilde{\omega})$')
# plt.title(fr'Dispersion relation: $\alpha={a}$, $\beta={b}$, $\Omega={om_c}$')
# plt.legend()

# # Add text box
# stability_text = "System: stable" if not instability_criterion else "System: unstable"
# plt.text(0.05, 0.95, stability_text, transform=plt.gca().transAxes, fontsize=12,
#          verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# plt.show()

# Define the complex function
def complex_function(x):
    return np.exp(1j * x)  # Example: e^(ix)

# Generate input values
x = np.linspace(0, 2 * np.pi, 500)  # Values from 0 to 2π

# Evaluate the function
y = complex_function(x)

# Extract real and imaginary parts
y_real = np.real(y)
y_imag = np.imag(y)

# Plot the real and imaginary parts
plt.figure(figsize=(8, 6))
plt.plot(x, y_real, label="Real Part", color="blue")
plt.plot(x, y_imag, label="Imaginary Part", color="orange")
plt.axhline(0, color="black", linewidth=0.8, linestyle="--")
plt.title("Real and Imaginary Parts of a Complex Function")
plt.xlabel("x")
plt.ylabel("Function Value")
plt.legend()
plt.grid(True)
plt.show()
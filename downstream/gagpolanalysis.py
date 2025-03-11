#data analysis for gag-pol ratios
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Generate hypothetical data points
n_points = 20
x = np.linspace(0.5, 2, n_points)  # gag/pol reads per kb ratios
y = 1 / x  # gag/pol adenine sequences per kb ratios (inverse proportion)

# Add some random noise to make it more realistic
np.random.seed(42)  # for reproducibility
x *= np.random.normal(1, 0.1, n_points)
y *= np.random.normal(1, 0.1, n_points)

# Calculate logarithms
log_x = np.log(x)
log_y = np.log(y)

# Calculate Pearson correlation coefficient
correlation_coefficient, p_value = stats.pearsonr(log_x, log_y)

# Create the plot
plt.figure(figsize=(12, 8))
plt.scatter(log_x, log_y, color='blue', label='Data points')

# Linear regression
coeffs = np.polyfit(log_x, log_y, 1)
poly = np.poly1d(coeffs)
plt.plot(log_x, poly(log_x), color='red', linestyle='--',
         label=f'Linear fit (slope: {coeffs[0]:.4f})')

# Labels and title
plt.xlabel('log(gag/pol reads per kb ratio)')
plt.ylabel('log(gag/pol adenine sequences per kb ratio)')
plt.title('Hypothetical Log-Log Plot of gag/pol Ratios')
plt.legend()
plt.grid(True, which="both", ls="-", alpha=0.2)

# Show the plot
plt.tight_layout()
plt.show()

# Print results
print(f"Correlation coefficient between log(x) and log(y): {correlation_coefficient:.6f}")
print(f"P-value: {p_value:.6f}")
print(f"Slope of the linear fit: {coeffs[0]:.4f}")
print(f"Intercept of the linear fit: {coeffs[1]:.4f}")

# Calculate R-squared
residuals = log_y - poly(log_x)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((log_y - np.mean(log_y))**2)
r_squared = 1 - (ss_res / ss_tot)
print(f"R-squared: {r_squared:.4f}")
import matplotlib.pyplot as plt
import numpy as np

# Data
data = {
    'A11_scPatho': {'A_seq_ratio': 18.2328190743338 / 12.4020, 'Reads_ratio': 2500 / 1715.28751753156},
    'A11_10x': {'A_seq_ratio': 18.2328190743338 / 12.4020, 'Reads_ratio': 4721.93211488251 / 3321.52875175316},
    'A10_scPatho': {'A_seq_ratio': 16.4796633941094 / 11.0966057441253, 'Reads_ratio': 2500 / 1715.28751753156},
    'A10_10x': {'A_seq_ratio': 16.4796633941094 / 11.0966057441253, 'Reads_ratio': 4721.93211488251 / 3321.52875175316}
}

# Create the plot
fig, ax = plt.subplots(figsize=(10, 10))

colors = {'A11': '#FF9999', 'A10': '#66B2FF'}
markers = {'scPatho': 'o', '10x': 's'}

for key, values in data.items():
    a_type, data_type = key.split('_')
    ax.scatter(values['A_seq_ratio'], values['Reads_ratio'],
               color=colors[a_type], marker=markers[data_type], s=100,
               label=f"{a_type} {data_type}")

# Plot the y=x line
x = np.linspace(1, 3, 100)
ax.plot(x, x, 'r--', alpha=0.5, label='y = x')

# Set labels and title
ax.set_xlabel('pol-A-sequences-per-kb/gag-A-sequences-per-kb')
ax.set_ylabel('gag-reads-per-kb/pol-reads-per-kb')
ax.set_title('Correlation between Adenine Sequences Ratio and Reads Ratio in HIV Gag-Pol')

# Customize the plot
ax.legend(loc='upper left')
ax.grid(True, linestyle='--', alpha=0.7)
ax.set_aspect('equal')
ax.set_xlim(1, 2)
ax.set_ylim(1, 2)

# Add text explaining the relationship
plt.text(0.05, 0.05,
         "The close alignment of data points with the y=x line demonstrates\n"
         "the 1-to-1 correlation between the pol/gag ratio of adenine\n"
         "sequences and the gag/pol ratio of reads. This strongly suggests\n"
         "that gag and pol are part of the same transcript, with differences\n"
         "in apparent expression levels due to transcript structure.",
         transform=ax.transAxes, fontsize=9, verticalalignment='bottom',
         bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.5))

plt.tight_layout()
plt.savefig('hiv_gag_pol_ratio_correlation_clean.png')
plt.close()

print("Clean plot has been saved as 'hiv_gag_pol_ratio_correlation_clean.png'")
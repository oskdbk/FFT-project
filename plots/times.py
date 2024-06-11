from tabulate import tabulate
import matplotlib.pyplot as plt
import numpy as np
lengths = [100, 1000, 10000, 25000, 50000, 100000]
data_1 = {
    "Sequential": [0.0006, 0.00824, 0.11576, 0.24376, 0.57128, 1.26604],
    "Inplace": [0.00036, 0.00604, 0.07588, 0.17732, 0.38432, 0.86768],
    "Parallel_1": [0.00004, 0.00668, 0.0808, 0.22068, 0.46476, 1.02504],
    "Parallel_2": [0.0042, 0.00492, 0.04916, 0.11136, 0.2502, 0.5514],
    "Parallel_4": [0.00152, 0.00372, 0.02944, 0.06152, 0.13532, 0.3084],
    "Parallel_8": [0.00304, 0.00388, 0.0176, 0.03744, 0.0774, 0.17268],
    "Parallel_12": [0.00388, 0.00484, 0.01836, 0.04492, 0.06028, 0.17736],
    "Parallel_16": [0.00508, 0.0054, 0.01844, 0.03708, 0.07416, 0.16936],
    "Parallel_32": [0.01052, 0.01064, 0.02008, 0.03904, 0.06516, 0.15632]
}




data_2 = {
    "Radix2": [0.0002, 0.00192, 0.04136, 0.082, 0.14492, 0.35928],
    "Radix2_Parallel_1": [0.00016, 0.00164, 0.04172, 0.07708, 0.14132, 0.3288],
    "Radix2_Parallel_2": [0.00044, 0.002, 0.02456, 0.04556, 0.09, 0.19292],
    "Radix2_Parallel_4": [0.00004, 0.0012, 0.01664, 0.02908, 0.05948, 0.12188],
    "Radix2_Parallel_8": [0.001, 0.00128, 0.01204, 0.02164, 0.0442, 0.08396],
    "Radix2_Parallel_12": [0.00124, 0.00176, 0.01124, 0.01932, 0.04104, 0.09156],
    "Radix2_Parallel_16": [0.00144, 0.0012, 0.01096, 0.01844, 0.0368, 0.0844],
    "Radix2_Parallel_32": [0.00292, 0.00356, 0.0106, 0.01732, 0.03364, 0.07536]
}


# Create a figure with a single subplot
fig, ax = plt.subplots()

colors = plt.get_cmap('viridis')(np.linspace(0, 1, len(data_2)))
for i, (l, temps) in enumerate(data_2.items()):
    color = colors[i] if i > 0 else "red"
    ax.plot(lengths, temps, label=l, color=color, linewidth=1)

# Set the labels and title
ax.set_xlabel('Length')
ax.set_ylabel('Time (in s)')
ax.set_title('Radix2 Time Comparison')
ax.set_xscale('log')
ax.set_yscale('log')

# Set log-log scale

ax.legend()

# Show the figure
plt.savefig("plot_radix2_loglog.png")

# Create a figure with a single subplot
fig, ax = plt.subplots()

colors = plt.get_cmap('viridis')(np.linspace(0, 1, len(data_1)))
for i, (l, temps) in enumerate(data_1.items()):
    color = colors[i] if i > 1 else ("red" if i>0 else "blue")
    ax.plot(lengths, temps, label=l, color=color, linewidth=1)

# Set the labels and title
ax.set_xlabel('Length')
ax.set_ylabel('Time (in s)')
ax.set_title('General Algorithm Time Comparison')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()

# Show the figure
plt.savefig("plot_general_loglog.png")

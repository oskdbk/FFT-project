from tabulate import tabulate
import numpy as np

data = [
    ["Length", "Sequential", "Inplace", "Parallel_1", "Parallel_2", "Parallel_4", "Parallel_8", "Parallel_12", "Parallel_16", "Parallel_32", "Natali"],
    [100, 0.0006, 0.0004, 0.0006, 0.0014, 0.0012, 0.0022, 0.0034, 0.0042, 0.0088, 0.0030],
    [1000, 0.0112, 0.0062, 0.0080, 0.0060, 0.0044, 0.0010, 0.0068, 0.0050, 0.0086, 0.0082],
    [10000, 0.1346, 0.0954, 0.1084, 0.0566, 0.0308, 0.0252, 0.0254, 0.0236, 0.0226, 0.0606],
    [25000, 0.3194, 0.2504, 0.2602, 0.1372, 0.0782, 0.0558, 0.0522, 0.0500, 0.0476, 0.2150],
    [50000, 0.7320, 0.5330, 0.5828, 0.3048, 0.1638, 0.1100, 0.0796, 0.0978, 0.0984, 0.5546],
    [100000, 1.6118, 1.2040, 1.3206, 0.6730, 0.3590, 0.2380, 0.2132, 0.2532, 0.2342, 1.6102]
]

import matplotlib.pyplot as plt

data2 = [
    ["Length", "Sequential", "Inplace", "2Radix", "Parallel_1", "Parallel_2", "Parallel_4", "Parallel_8", "Parallel_12", "Parallel_16", "Parallel_32", "Natali"],
    [100,     0.00045, 0.00025, 0.00015, 0.0003,  0.001,   0.00125, 0.00225, 0.00345, 0.00435, 0.0096,  0.00365],
    [1000,    0.0066,  0.0045,  0.0015,  0.00505, 0.004,   0.0031,  0.0033,  0.00405, 0.0051,  0.0101,  0.0093],
    [10000,   0.08385, 0.0571,  0.03035, 0.06245, 0.03365, 0.01915, 0.0142,  0.0145,  0.0148,  0.01595, 0.04585],
    [25000,   0.1932,  0.1424,  0.06225, 0.1517,  0.09095, 0.05175, 0.03085, 0.02725, 0.03005, 0.0301,  0.13915],
    [50000,   0.4534,  0.31985, 0.12875, 0.35855, 0.2057,  0.1264,  0.07735, 0.07995, 0.07765, 0.08005, 0.3952],
    [100000,  0.912,   0.653,   0.2625,  0.7875,  0.4572,  0.2902,  0.1762,  0.1766,  0.1768,  0.177,   0.8545]
]

lengths = [row[0] for row in data2[1:]]
sequential = [row[1] for row in data2[1:]]
inplace = [row[2] for row in data2[1:]]
radix = [row[3] for row in data2[1:]]
parallel_1 = [row[4] for row in data2[1:]]
parallel_2 = [row[5] for row in data2[1:]]
parallel_4 = [row[6] for row in data2[1:]]
parallel_8 = [row[7] for row in data2[1:]]
parallel_12 = [row[8] for row in data2[1:]]
parallel_16 = [row[9] for row in data2[1:]]
parallel_32 = [row[10] for row in data2[1:]]
natali = [row[11] for row in data2[1:]]

plt.plot(lengths, sequential, label="Sequential", color="darkred")
#plt.plot(lengths, inplace, label="Inplace", color = "red")
plt.plot(lengths, radix, label="2Radix", color = "orange")
colors = np.linspace(0, 1, 7)
colormap = plt.cm.viridis

#plt.plot(lengths, parallel_1, label="Parallel_1", color=colormap(colors[0]))
plt.plot(lengths, parallel_2, label="Parallel_2", color=colormap(colors[1]))
plt.plot(lengths, parallel_4, label="Parallel_4", color=colormap(colors[2]))
#plt.plot(lengths, parallel_8, label="Parallel_8", color=colormap(colors[3]))
plt.plot(lengths, parallel_12, label="Parallel_12", color=colormap(colors[4]))
#plt.plot(lengths, parallel_16, label="Parallel_16", color=colormap(colors[5]))
plt.plot(lengths, parallel_32, label="Parallel_32", color=colormap(colors[6]))
#plt.plot(lengths, natali, label="Natali")

plt.xlabel("Length")
plt.ylabel("Time (s)")
#plt.xscale("log")
#plt.yscale("log")
plt.title("Execution Time vs Length")
plt.legend()
plt.grid(True)
plt.show()

print(tabulate(data2[1:], headers=data[0], floatfmt=".4f", tablefmt="grid"))
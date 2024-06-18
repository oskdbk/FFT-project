from tabulate import tabulate
import matplotlib.pyplot as plt
import numpy as np
lengths = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000, 23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 31000, 32000, 33000, 34000, 35000, 36000, 37000, 38000, 39000, 40000]
until = 28
# read data from file
with open("timing/general_times.csv", "r") as f:
    headers = f.readline().strip().split(",")[:-1]
    data_1 = {}
    for h in headers:
        data_1[h] = []
    for line in f:
        line = line.strip().split(",")[:-1]
        for i, val in enumerate(line):
            data_1[headers[i]].append(float(val))

# do the same for radix2
with open("timing/radix_times.csv", "r") as f:
    headers = f.readline().strip().split(",")[:-1]
    data_2 = {}
    for h in headers:
        data_2[h] = []
    for line in f:
        line = line.strip().split(",")[:-1]
        for i, val in enumerate(line):
            data_2[headers[i]].append(float(val))

# print(data_1)
# print(data_2)

# data_1 = {
#     "Sequential": [0.0006, 0.00824, 0.11576, 0.24376, 0.57128, 1.26604],
#     "Inplace": [0.00036, 0.00604, 0.07588, 0.17732, 0.38432, 0.86768],
#     "Parallel_1": [0.0004, 0.00668, 0.0808, 0.22068, 0.46476, 1.02504],
#     "Parallel_2": [0.0042, 0.00492, 0.04916, 0.11136, 0.2502, 0.5514],
#     "Parallel_4": [0.00152, 0.00372, 0.02944, 0.06152, 0.13532, 0.3084],
#     "Parallel_8": [0.00304, 0.00388, 0.0176, 0.03744, 0.0774, 0.17268],
#     "Parallel_12": [0.00388, 0.00484, 0.01836, 0.04492, 0.06028, 0.17736],
#     "Parallel_16": [0.00508, 0.0054, 0.01844, 0.03708, 0.07416, 0.16936],
#     "Parallel_32": [0.01052, 0.01064, 0.02008, 0.03904, 0.06516, 0.15632]
# }




# data_2 = {
#     "Radix2": [0.0002, 0.00192, 0.04136, 0.082, 0.14492, 0.35928],
#     "Radix2_Parallel_1": [0.00016, 0.00164, 0.04172, 0.07708, 0.14132, 0.3288],
#     "Radix2_Parallel_2": [0.00044, 0.002, 0.02456, 0.04556, 0.09, 0.19292],
#     "Radix2_Parallel_4": [0.0004, 0.0012, 0.01664, 0.02908, 0.05948, 0.12188],
#     "Radix2_Parallel_8": [0.001, 0.00128, 0.01204, 0.02164, 0.0442, 0.08396],
#     "Radix2_Parallel_12": [0.00124, 0.00176, 0.01124, 0.01932, 0.04104, 0.09156],
#     "Radix2_Parallel_16": [0.00144, 0.0012, 0.01096, 0.01844, 0.0368, 0.0844],
#     "Radix2_Parallel_32": [0.00292, 0.00356, 0.0106, 0.01732, 0.03364, 0.07536]
# }

for until in [10,19,28]:
    break

    # Create a figure with a single subplot
    fig, ax = plt.subplots()

    colors = plt.get_cmap('viridis')(np.linspace(0, 1, len(data_2)))
    for i, (l, temps) in enumerate(data_2.items()):
        color = colors[i] if i > 0 else "red"
        ax.plot(lengths[:until], temps[:until], label=l, color=color, linewidth=1)

    # Set the labels and title
    ax.set_xlabel('Length')
    ax.set_ylabel('Time (in s)')
    ax.set_title('Radix2 Time Comparison')
    # ax.set_xscale('log')
    # ax.set_yscale('log')

    # Set log-log scale

    ax.legend()

    # Show the figure
    plt.savefig(f"plot_radix2_until_{until}.png")
    #plt.show()

    # Create a figure with a single subplot
    fig, ax = plt.subplots()

    colors = plt.get_cmap('viridis')(np.linspace(0, 1, len(data_1)))
    for i, (l, temps) in enumerate(data_1.items()):
        color = colors[i] if i > 1 else ("red" if i>0 else "blue")
        ax.plot(lengths[:until], temps[:until], label=l, color=color, linewidth=1)

    # Set the labels and title
    ax.set_xlabel('Length')
    ax.set_ylabel('Time (in s)')
    ax.set_title('General Algorithm Time Comparison')
    # ax.set_xscale('log')
    # ax.set_yscale('log')


    ax.legend()

    # Show the figure
    plt.savefig(f"plot_general_until_{until}.png")


# plot radix vs general
for until in [10,19,28, -1]:
    continue
    fig, ax = plt.subplots()
    general_colors = plt.get_cmap('Greens')(np.linspace(0.5, 0.8, 2))
    radix_colors = plt.get_cmap('Oranges')(np.linspace(0.5, 0.8, 2))
    c = plt.get_cmap('viridis')
    temps = data_1["Inplace"]
    ax.plot(lengths[:until], temps[:until], label="general inplace", color = c(0.2), linewidth=1, linestyle="dashed")
    temps = data_2["Radix2"]
    ax.plot(lengths[:until], temps[:until], label="radix serial", color = c(0.8), linewidth=1, linestyle="dashed")
    temps = data_1["Parallel_12"]
    ax.plot(lengths[:until], temps[:until], label="general parallel 12", color = c(0.2), linewidth=1)
    temps = data_2["Radix2_Parallel_12"]
    ax.plot(lengths[:until], temps[:until], label="radix parallel 32", color = c(0.8), linewidth=1)



    # Set the labels and title
    ax.set_xlabel('Length')
    ax.set_ylabel('Time (in s)')
    ax.set_title('Comparison between the two algorithms')

    ax.legend()
    plt.savefig(f"plot_comparison_until_{until}.png")


# # create a table
# all_data = {
#     "radix": {
#         "sequential": data_2["Radix2"],
#         "parallel_2": data_2["Radix2_Parallel_2"],
#         "parallel_4": data_2["Radix2_Parallel_4"],
#         #"parallel_8": data_2["Radix2_Parallel_8"],
#         "parallel_12": data_2["Radix2_Parallel_12"],
#         #"parallel_16": data_2["Radix2_Parallel_16"],
#         "parallel_32": data_2["Radix2_Parallel_32"]
#     },
#     "general": {
#         "sequential": data_1["Sequential"],
#         "inplace": data_1["Inplace"],
#         "parallel_2": data_1["Parallel_2"],
#         "parallel_4": data_1["Parallel_4"],
#         #"parallel_8": data_1["Parallel_8"],
#         "parallel_12": data_1["Parallel_12"],
#         #"parallel_16": data_1["Parallel_16"],
#         "parallel_32": data_1["Parallel_32"]
#     }
# }

# # create a table
# table = []
# for i, length in enumerate(lengths):
#     row = [length]
#     for alg in all_data:
#         for parallel in all_data[alg]:
#             val = all_data[alg][parallel][i]
#             row.append(f"{val:.5f}")
#             #row.append(all_data[alg][parallel][i])
#     table.append(row)
#     row2 = [" "]
#     for alg in all_data:
#         best = min([all_data[alg][parallel][i] for parallel in all_data[alg]])
#         worst = max([all_data[alg][parallel][i] for parallel in all_data[alg]])
#         for parallel in all_data[alg]:
#             val = all_data[alg][parallel][i]
#             try:
#                 ref = all_data[alg]["inplace"][i]
#             except:
#                 ref = all_data[alg]["sequential"][i]
#             ratio = val/ref
#             if ratio == 1:
#                 row2.append(f" ")
#             elif val == worst:
#                 row2.append(f"\\textcolor{{worst}}{{{ratio:.2f}}}")
#             elif val == best:
#                 row2.append(f"\\textcolor{{best}}{{{ratio:.2f}}}")
#             else:
#                 row2.append(f"\\textcolor{{greyish}}{{{ratio:.2f}}}")
#     table.append(row2)
# # make into latex table
# l = ""
# for row in table:
#     l += " & ".join([str(i) for i in row]) + " \\\\ \n"

# print(l)

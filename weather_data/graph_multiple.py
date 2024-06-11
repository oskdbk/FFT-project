import numpy as np

import matplotlib.pyplot as plt

# Read the data from the file
temperatures = []
with open('weather_data.csv', 'r') as file:
    for line in file:
        temperature = float(line.strip())
        temperatures.append(temperature)

# Read the data from the cleaned file

with open('weather_cleaned_radix.csv', 'r') as file:
    labels = file.readline().split(",")[:-1]
    cleaned_temperatures = {l: [] for l in labels}
    for line in file:
        temps = line.split(",")[:-1]
        for i,l in enumerate(labels):
            temperature = float(temps[i])
            cleaned_temperatures[l].append(temperature)

# Create a figure with a single subplot
fig, ax = plt.subplots()

# Plot the temperatures and cleaned temperatures on the same subplot
ax.scatter(range(len(temperatures)), temperatures, label='temperature data', color='lightblue', s=6)
# Generate a gradient of colors
colors = plt.get_cmap('viridis')(np.linspace(0, 1, len(cleaned_temperatures)))
for i, (l, temps) in enumerate(cleaned_temperatures.items()):
    ax.plot(range(len(temps)), temps, label=l+" terms", color=colors[i], linewidth=1)

# Set the labels and title
ax.set_xlabel('Time')
ax.set_ylabel('Temperature')
ax.set_title('Temperature Data Reconstruction with 2Radix')

# Add a legend
ax.legend()

# Show the figure
plt.savefig("plot_multiple_2radix.png")
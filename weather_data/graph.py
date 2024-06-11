import matplotlib.pyplot as plt

# Read the data from the file
temperatures = []
with open('weather_data.csv', 'r') as file:
    for line in file:
        temperature = float(line.strip())
        temperatures.append(temperature)

# Read the data from the cleaned file
cleaned_temperatures = []
with open('weather_data_cleaned.csv', 'r') as file:
    for line in file:
        temperature = float(line.strip())
        cleaned_temperatures.append(temperature)

# Create a figure with a single subplot
fig, ax = plt.subplots()

# Plot the temperatures and cleaned temperatures on the same subplot
ax.scatter(range(len(temperatures)), temperatures, label='Temperature Data', color='lightblue', s=6)
ax.plot(cleaned_temperatures, label='Cleaned Temperature Data', color='lightcoral', linewidth=2)

# Set the labels and title
ax.set_xlabel('Time')
ax.set_ylabel('Temperature')
ax.set_title('Temperature Data')

# Add a legend
ax.legend()

# Show the figure
plt.show()

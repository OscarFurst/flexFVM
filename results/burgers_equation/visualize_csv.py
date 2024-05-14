import os
import glob

import matplotlib.pyplot as plt
import matplotlib.animation as animation

import numpy as np

# Path to the directory containing the CSV files
csv_dir = './burger_time_series'

# Get a list of all CSV files in the directory
csv_files = glob.glob(os.path.join(csv_dir, 'time*.csv'))

# Sort the CSV files based on their numerical suffix
csv_files.sort(key=lambda x: float(
    os.path.splitext(os.path.basename(x))[0][4:]))

# Get the array of time data from the filenames
times = np.array([float(os.path.splitext(os.path.basename(x))[0][4:])
                  for x in csv_files])

# Create a figure and axis for the plot
fig, ax = plt.subplots()
ax.set_title('1D Burgers Equation with periodic boundary conditions')
ax.set_xlabel('x / m')
ax.set_ylabel('u / m.s-1')

# Load the initial data to set the axis limits
data = np.loadtxt(csv_files[0], delimiter=',', skiprows=1)
data = data[data[:, 0].argsort()]
ax.set_xlim(data[:, 0].min(), data[:, 0].max())
ax.set_ylim(data[:, 3].min(), data[:, 3].max())

# Initialize an empty line object for the plot
line, = ax.plot([], [], lw=2)
timelabel = ax.text(-6, 0.9, '', fontsize=15, ha='left', va='top')

# Define the initialization function for the animation


def init():
    line.set_data([], [])
    timelabel.set_text('')
    return line, timelabel

# Define the update function for the animation


def update(frame):
    # Load the data from the current CSV file
    data = np.loadtxt(csv_files[frame], delimiter=',', skiprows=1)

    # Update the line object with the new data
    data = data[data[:, 0].argsort()]
    line.set_data(data[:, 0], data[:, 3])
    timelabel.set_text('Time: {:.2f} s'.format(times[frame]))
    return line, timelabel


# Create the animation
ani = animation.FuncAnimation(fig, update, frames=len(
    csv_files), init_func=init, blit=True, interval=200/5)

# Show the animation
plt.show()

# Save the animation as an MP4 file
ani.save('burger_time_series.mp4', writer='ffmpeg')

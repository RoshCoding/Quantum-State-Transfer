import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

# Load data from qst-out-n.txt in the same folder as this script
script_dir = os.path.dirname(__file__)
data_path = os.path.join(script_dir, 'qst-out-p.txt')
print(f"Loading data from: {data_path}")
data = np.loadtxt(data_path)

num_frames, num_points = data.shape
x = np.arange(num_points)

fig, ax = plt.subplots()
line, = ax.plot(x, data[0])
ax.set_ylim(np.min(data), np.max(data))
ax.set_xlim(0, num_points - 1)
ax.set_xlabel('Position')
ax.set_ylabel('Probability')
ax.set_title('Quantum State Evolution')

def update(frame):
    line.set_ydata(data[frame])
    return line,

ani = animation.FuncAnimation(fig, update, frames=num_frames, blit=True, interval=50)
output_path = os.path.join(script_dir, 'out-p.mp4')
ani.save(output_path, writer='ffmpeg', fps=10)
plt.close(fig)
print('Animation saved as out-p.mp4')
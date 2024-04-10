import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'font.size': 14})

opts = np.array([0, 1, 2, 3])

data = np.array([
    # I19 Serial
    [1.27575, 0.19651, 0.05377, 0.05372],
    # I19 Parallel
    [0.12185, 0.02686, 0.01913, 0.01516],
    # I24 Serial
    [0.98183, 0.19668, 0.13139, 0.13204],
    # I24 Parallel
    [0.10101, 0.02678, 0.01868, 0.01958],
    # G11 Serial
    [1.00503, 0.29993, 0.29457, 0.07263],
    # G11 Parallel
    [0.06986, 0.02349, 0.02156, 0.01541]
])

fig, axs = plt.subplots(1, 2)

axs[0].plot(opts, data[0], '.-r', label='Intel 19', ms=16, lw=3)
axs[0].plot(opts, data[2], '.-g', label='Intel 24', ms=16, lw=3)
axs[0].plot(opts, data[4], '.-b', label='GNU 11', ms=16, lw=3)

axs[1].plot(opts, data[1], '.-r', label='Intel 19', ms=16, lw=3)
axs[1].plot(opts, data[3], '.-g', label='Intel 24', ms=16, lw=3)
axs[1].plot(opts, data[5], '.-b', label='GNU 11', ms=16, lw=3)

axs[0].set_ylabel('Runtime (s)')
axs[0].set_xlabel('Optimization Level')
axs[0].set_xticks(opts)
axs[0].set_title('Serial')
axs[1].set_xlabel('Optimization Level')
axs[1].set_xticks(opts)
axs[1].set_title('Parallel')
axs[1].legend()

plt.show()
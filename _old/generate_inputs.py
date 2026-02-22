import scipy.ndimage
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from util import *

np.random.seed(1)


############

# m = 20
# n = 10
#
# N = 10
# a = 3.0
# x = np.arange(1, N+1)
# weights = x ** (-a)
# weights /= weights.sum()
# bounded_zipf = stats.rv_discrete(name='bounded_zipf', values=(x, weights))
# sample = bounded_zipf.rvs(size=(m * n))
# rewards = np.reshape(sample, (m, n)) - 1.0
#
# plt.imshow(rewards, cmap='jet', interpolation='bilinear')
# plt.savefig('plot_heatmap_test.pdf')

rewards = read_input_from_file_id(1017)
rewards = scipy.ndimage.zoom(rewards, 0.33, order=1)
[m, n] = rewards.shape
plt.imshow(rewards, cmap='jet', interpolation='bilinear')
plt.savefig('plot_heatmap_test.pdf')


m = 5
n = 50
out = np.zeros((m, n))

out[0][1] = 3000
out[0][n-1] = 3000

out[1][1] = 300
out[1][n-1] = 300

out[2][1] = 20
out[2][n-1] = 20

out[3][1] = 300
out[3][n-1] = 300

out[4][1] = 3000
out[4][n-1] = 3000

np.savetxt("input/input-0.csv", out, delimiter=",", fmt='%.2f')
azz = 0

############


a_vec = [0.0, 0.0, 0.8, 0.8, 1.9, 1.9, 2.7, 2.7]
size_vec = [[100, 50], [50, 100], [100, 50], [50, 100], [100, 50], [50, 100], [100, 50], [50, 100]]
base_vec = [100, 200, 300, 400, 500, 600, 700, 800]

for i in range(0, len(base_vec)):
    for j in range(0, 30):
        input = base_vec[i] + j

        m = size_vec[i][0]
        n = size_vec[i][1]

        N = 100
        a = a_vec[i]
        x = np.arange(1, N+1)
        weights = x ** (-a)
        weights /= weights.sum()
        bounded_zipf = stats.rv_discrete(name='bounded_zipf', values=(x, weights))

        scale = 5
        sample = bounded_zipf.rvs(size=(m/scale * n/scale))
        rewards = np.reshape(sample, (m/scale, n/scale))
        rewards = scipy.ndimage.zoom(rewards, scale, order=1) - 1.0
        rewards = np.array(rewards, dtype=float)

        # append zeros first-last columns
        tmp = np.copy(rewards)
        rewards = np.zeros((m, n+2))
        for k in range(1, n+1):
            rewards[:, k] = tmp[:, k-1]

        np.savetxt("input/input-%d.csv" % input, rewards, delimiter=",", fmt='%.2f')

        if j == 0:
            plt.imshow(rewards, cmap='jet', interpolation='bilinear')
            plt.savefig('plot_heatmap_%d.pdf' % base_vec[i])

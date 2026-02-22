import csv
import numpy as np
from scipy.stats import sem, t
from scipy import mean


def read_input_from_file_id(input):
    reader = csv.reader(open("input/input-%d.csv" % input, "r"), delimiter=",")
    x = list(reader)
    rewards = np.array(x).astype("float")

    return rewards


def get_confidence_interval(data, confidence):
    size = np.size(data, 1)
    m = np.zeros(size)
    h = np.zeros(size)

    for i in range(0, size):
        current_data = data[:, i]
        n = len(current_data)
        m[i] = mean(current_data)
        std_err = sem(current_data)
        h[i] = std_err * t.ppf((1 + confidence) / 2, n - 1)

    return [m, h]

# input = 1001
# R = read_input_from_file_id(input)
# [m, n] = R.shape
# tmp = np.copy(R)
# rewards = np.zeros((m, n+2))
# for j in range(1, n+1):
#     rewards[:, j] = tmp[:, j-1]
#
# np.savetxt("input/input-%d.csv" % input, rewards, delimiter=",", fmt='%.2f')


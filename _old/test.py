from algorithms import *
import matplotlib.pyplot as plt


# for input in range(1010, 1020):
#     rewards = read_input_from_file_id(input)
#     plt.imshow(rewards, cmap='jet', interpolation='bilinear')
#     plt.savefig('plot_heatmap_%d.pdf' % input)
#

input = 0
budget = 120
rewards = read_input_from_file_id(input)
#
# m = 30
# n = 100
# rewards = np.zeros((m, n+2))
# for i in range(0, m-1):
#     rewards[i][1] = 100*(i+2)-1
# rewards[m-1][m-2] = 100*((m-1) + 2*(m-2))
#
# np.savetxt("input/input-%d.csv" % input, rewards, delimiter=",", fmt='%.2f')

# q = []
# base = 600
# for input in range(base, base+30):
#     rewards = read_input_from_file_id(input)
#     azz = rewards[rewards>5]
#     azz = 100*len(azz)/(50*100.0)
#     q.append(azz)
# print(np.mean(q))




# out = greedy_full_row(input, budget)
# print("greedy_full_row (gfr) -> reward=%.2f, cost=%d/%d" % (out["reward"], out["cost"], budget))
#
out = opt_full_row(input, budget)
print("opt_full_row (ofr) -> reward=%.2f, cost=%d/%d" % (out["reward"], out["cost"], budget))

# out = opt_partial_row_single_column_cpp(input, budget)
# print("opt_partial_row_single_column (oprsc) -> reward=%.2f, cost=%d/%d" % (out["reward"], out["cost"], budget))
#
out = greedy_partial_row(input, budget)
print("greedy_partial_row (gpr) -> reward=%.2f, cost=%d/%d" % (out["reward"], out["cost"], budget))

out = heuristic_partial_row(input, budget)
print("heuristic_partial_row (hprgc) -> reward=%.2f, cost=%d/%d" % (out["reward"], out["cost"], budget))





# out = heuristic_1(input, budget)
# print("heuristic_1 -> reward=%.2f, cost=%d/%d" % (out["reward"], out["cost"], budget))
# out = heuristic_2(input, budget)
# print("heuristic_2 -> reward=%.2f, cost=%d/%d" % (out["reward"], out["cost"], budget))
# out = heuristic_3(input, budget)
# print("heuristic_3 -> reward=%.2f, cost=%d/%d" % (out["reward"], out["cost"], budget))
# out = heuristic_4(input, budget)
# print("heuristic_4 -> reward=%.2f, cost=%d/%d" % (out["reward"], out["cost"], budget))
# out = heuristic_5(input, budget)
# print("heuristic_5 -> reward=%.2f, cost=%d/%d" % (out["reward"], out["cost"], budget))





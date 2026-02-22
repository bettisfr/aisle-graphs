from algorithms import *
import matplotlib.pyplot as plt
from util import *


colors = ['dodgerblue', 'blue', 'black', 'lightcoral', 'red']
styles = ['dashed', 'solid', 'solid', 'dashed', 'solid']
algs = ['gfr', 'ofr', 'oprsc', 'gpr', 'hprgc']
algs_long = ['GFr', 'OFr', 'OSc', 'GPr', 'HGc']
alg_plot = [True, True, True, True, True]

base = 1010
exps = 10

[m, n] = read_input_from_file_id(base).shape
max_budget = n * m + 3 * m
budgets = np.linspace(0, max_budget, num=20)
size = len(budgets)

rewards = np.zeros((len(algs), exps, len(budgets)))
costs = np.zeros((len(algs), exps, len(budgets)))

for j in range(0, exps):
    input = base + j

    for k in range(0, len(algs)):
        if not alg_plot[k]:
            continue

        alg = algs[k]

        file_path = 'output/%s/out_%s-%d.txt' % (alg, alg, input)
        with open(file_path) as fp:
            for cnt, line in enumerate(fp):
                if cnt < len(budgets):
                    data = line.strip().split(",")
                    costs[k][j][cnt] = float(data[1])
                    rewards[k][j][cnt] = float(data[2])

    max_rewards = np.sum(read_input_from_file_id(input))
    rewards[:,j,:] = 100 * rewards[:,j,:] / max_rewards

budgets = 100 * np.array(budgets) / budgets[-1]
np.savetxt("data_b.csv", budgets, delimiter=",", fmt='%.2f')

plt.clf()
fig = plt.figure(figsize=(2.5, 2.5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xlabel('$B (\\%)$')
plt.ylabel('$R (\\%)$')
plt.title('$G = [%d \\times %s]$' % (m, n-2))

axes = plt.gca()
axes.set_xlim([0, 100])
axes.set_ylim([0, 100])

for k in range(0, len(algs)):
    if not alg_plot[k]:
        continue

    [avg, std] = get_confidence_interval(rewards[k], 0.95)
    np.savetxt("data_%d_avg.csv" % (k), avg, delimiter=",", fmt='%.2f')
    np.savetxt("data_%d_std.csv" % (k), std, delimiter=",", fmt='%.2f')
    plt.errorbar(budgets, avg, std, c=colors[k], label=algs_long[k], linewidth=0.75, linestyle=styles[k])

plt.legend(loc='lower right', prop={'size': 6})

plt.tight_layout(pad=0.5)
plt.savefig('plot_algs_comparison_real_%d.pdf' % base)

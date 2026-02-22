from algorithms import *
import matplotlib.pyplot as plt
from util import *


colors = ['dodgerblue', 'blue', 'black', 'lightcoral', 'red']
styles = ['dashed', 'solid', 'solid', 'dashed', 'solid']
algs = ['gfr', 'ofr', 'oprsc', 'gpr', 'hprgc']
algs_long = ['GFr', 'OFrI', 'OptSc', 'GPr', 'HGc']
alg_plot = [True, True, True, True, True]

base_vec = [100, 200, 300, 400, 500, 600, 700, 800]
a_vec = [0.0, 0.0, 0.8, 0.8, 1.9, 1.9, 2.7, 2.7]

base_vec = [1010]
a_vec = [0.0]

n_exps = 10

for i in range(0, len(base_vec)):
    [m, n] = read_input_from_file_id(base_vec[i]).shape
    max_budget = n * m + 3 * m
    budgets = np.linspace(0, max_budget, num=20)
    size = len(budgets)

    rewards = np.zeros((len(algs), n_exps, len(budgets)))
    costs = np.zeros((len(algs), n_exps, len(budgets)))

    for j in range(0, n_exps):
        input = base_vec[i] + j

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
    # np.savetxt("data_%d_b.csv" % i, budgets, delimiter=",", fmt='%.2f')

    plt.clf()
    fig = plt.figure(figsize=(2.1, 2.1))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel('$B (\\%)$')
    plt.ylabel('$R (\\%)$')
    plt.title('$A = [%d \\times %s]$' % (m, n-2))

    axes = plt.gca()
    axes.set_xlim([0, 100])
    axes.set_ylim([0, 100])

    axes = plt.gca()

    major_ticks = np.arange(0, 100 + 1, 20)
    minor_ticks = np.arange(0, 100 + 1, 5)

    axes.set_xticks(major_ticks)
    axes.set_xticks(minor_ticks, minor=True)
    axes.set_yticks(major_ticks)
    axes.set_yticks(minor_ticks, minor=True)

    # And a corresponding grid
    axes.grid(which='both')

    # Or if you want different settings for the grids:
    axes.grid(which='minor', alpha=0.1)
    axes.grid(which='major', alpha=0.25)

    axes.xaxis.set_label_coords(0.5, -0.17)
    axes.yaxis.set_label_coords(-0.17, 0.5)

    for k in range(0, len(algs)):
        if not alg_plot[k]:
            continue

        [avg, std] = get_confidence_interval(rewards[k], 0.95)
        # np.savetxt("data_%d_%d_avg.csv" % (i, k), avg, delimiter=",", fmt='%.2f')
        # np.savetxt("data_%d_%d_std.csv" % (i, k), std, delimiter=",", fmt='%.2f')
        plt.errorbar(budgets, avg, std, c=colors[k], label=algs_long[k], linewidth=0.75, linestyle=styles[k])

    plt.legend(loc='lower right', prop={'size': 6})

    plt.tight_layout(pad=0.2)
    plt.savefig('plot_algs_comparison_%d.pdf' % base_vec[i])

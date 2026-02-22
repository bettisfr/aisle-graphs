from algorithms import *
import matplotlib.pyplot as plt
from util import *


colors = ['purple', 'teal', 'darkorange', 'olive', 'teal']
styles = ['solid', 'solid', 'solid', 'solid', 'solid']
algs = ['gfr', 'ofr', 'oprsc', 'gpr', 'hprgc']
algs_long = ['GFr', 'OFr', 'OSc', 'GPr', 'HGc']

# base_vec = [100, 300, 500, 700]
base_vec = [200, 400, 600, 800]
a_vec = [0.0, 0.8, 1.9, 2.7]

for i in range(0, len(algs)):
    alg = algs[i]

    [m, n] = read_input_from_file_id(base_vec[0]).shape
    max_budget = n * m + 3 * m
    budgets = np.linspace(0, max_budget, num=20)
    size = len(budgets)

    rewards = np.zeros((len(base_vec), 30, len(budgets)))
    costs = np.zeros((len(base_vec), 30, len(budgets)))

    for j in range(0, len(base_vec)):
        for k in range(0, 30):
            input = base_vec[j] + k

            file_path = 'output/%s/out_%s-%d.txt' % (alg, alg, input)
            with open(file_path) as fp:
                for cnt, line in enumerate(fp):
                    if cnt < len(budgets):
                        data = line.strip().split(",")
                        costs[j][k][cnt] = float(data[1])
                        rewards[j][k][cnt] = float(data[2])

            max_rewards = np.sum(read_input_from_file_id(input))
            rewards[j,k,:] = 100 * rewards[j,k,:] / max_rewards

    budgets = 100 * np.array(budgets) / budgets[-1]

    plt.clf()
    fig = plt.figure(figsize=(2.5, 2.5))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel('$B (\\%)$')
    plt.ylabel('$R (\\%)$')
    plt.title('$G = [%d \\times %s], %s$' % (m, n-2, algs_long[i]))

    axes = plt.gca()
    axes.set_xlim([0, 100])
    axes.set_ylim([0, 100])

    for k in range(0, len(base_vec)):
        [avg, std] = get_confidence_interval(rewards[k], 0.95)

        azz = np.block([[budgets], [avg], [std]])
        azz = np.transpose(azz)
        np.savetxt("data_%d_%s.csv" % (i, k), azz, delimiter="   ", fmt='%.2f')

        plt.errorbar(budgets, avg, std, c=colors[k], label=a_vec[k], linewidth=0.75, linestyle=styles[k])

    plt.legend(loc='lower right', prop={'size': 6})

    plt.tight_layout(pad=0.5)
    plt.savefig('plot_algs_comparison_alg_%s.pdf' % algs[i])

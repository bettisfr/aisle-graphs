from algorithms import *


base_vec = [100, 200, 300, 400, 500, 600, 700, 800]

for i in range(0, len(base_vec)):
    for j in range(0, 30):
        input = base_vec[i] + j

        print("input=%d" % input)

        # # Our algorithms
        # print("#### opt_full_row (ofr)")
        # execute_batch('ofr', input)
        #
        # print("#### opt_partial_row_single_column (oprsc)")
        # execute_batch('oprsc', input)

        print("#### heuristic_partial_row (hprgc)")
        execute_batch('hprgc', input)

        # # Their algorithms
        # print("#### greedy_full_row (gfr)")
        # execute_batch('gfr', input)
        #
        # print("#### greedy_partial_row (gpr)")
        # execute_batch('gpr', input)



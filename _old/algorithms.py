import numpy as np
from math import *
import csv
import win32com.client
import subprocess
from util import *
import os


def heuristic_partial_row(input, budget):
    rewards = np.zeros(5)
    costs = np.zeros(5)

    out = heuristic_1(input, budget)
    rewards[0] = out["reward"]
    costs[0] = out["cost"]

    out = heuristic_2(input, budget)
    rewards[1] = out["reward"]
    costs[1] = out["cost"]

    out = heuristic_3(input, budget)
    rewards[2] = out["reward"]
    costs[2] = out["cost"]

    out = heuristic_4(input, budget)
    rewards[3] = out["reward"]
    costs[3] = out["cost"]

    out = heuristic_5(input, budget)
    rewards[4] = out["reward"]
    costs[4] = out["cost"]

    reward = np.max(rewards)
    cost = costs[np.argmax(rewards)]

    output = {"cycle": [], "reward": reward, "cost": cost}
    return output


# ofr(B) + osc(B_res)
def heuristic_1(input, budget):
    rewards = read_input_from_file_id(input)
    [m, n] = rewards.shape

    # hprgc algorithm
    output = opt_full_row(input, budget)
    reward_1 = output["reward"]
    cost_1 = output["cost"]

    S = output["traversed"]
    traversed = len(S)

    if traversed == 0:
        return output

    last_row = max(S)

    budget_res = budget - cost_1
    if budget_res % 2 == 1:
        budget_res = budget_res - 1

    if budget_res > 0:
        S_t = np.sort(S)

        # G_L
        tmp = np.copy(rewards[:, 0:(n/2)])
        G_L = np.zeros((last_row-traversed+1, np.size(tmp, 1)))
        tmp_i = 0
        for i in range(0, last_row):
            if i not in S_t:
                G_L[tmp_i] = tmp[i]
                tmp_i = tmp_i+1

        # G_L "lower left corner"
        needed_B = (m-last_row-1)*2
        cost_2 = budget_res
        if budget_res >= needed_B:
            # it is not free for "go and back" there...
            budget_res = budget_res - needed_B
            G_L_plus = np.zeros((m-last_row-1, np.size(tmp, 1)))
            tmp_i = 0
            for i in range(last_row+1, m):
                G_L_plus[tmp_i] = tmp[i]
                tmp_i = tmp_i + 1

            G_L = np.block([[G_L], [G_L_plus]])

        # G_R
        tmp = np.copy(rewards[:, (n/2):n])
        G_R = []
        for i in range(0, traversed-1, 2):
            for k in range(S_t[i]+1, S_t[i+1]):
                G_R.append(tmp[k])
        G_R = np.array(G_R)

        if len(G_R > 0):
            G_R = np.flip(G_R, 1)
            G = np.block([[G_L], [G_R]])
        else:
            G = np.copy(G_L)

        if len(G > 0):
            [m_G, n_G] = G.shape

            tmp_file = "input/tmp-G.csv"
            np.savetxt(tmp_file, G, delimiter=",", fmt='%.0f')

            fast = True
            if fast:
                output_sc = opt_partial_row_single_column_novertical_cpp(m_G, n_G, budget_res)
            else:
                output_sc = opt_partial_row_single_column_novertical(m_G, n_G, budget_res)

            reward_2 = output_sc["reward"]
            # cost_2 = output_sc["cost"]

            cost_tot = cost_1 + cost_2
            reward_tot = reward_1 + reward_2

            output["reward"] = reward_tot
            output["cost"] = cost_tot

    return output


# osc_L + osc_R => full + osc
def heuristic_2(input, budget):
    rewards = read_input_from_file_id(input)
    [m, n] = rewards.shape
    output = {"cycle": [], "reward": 0, "cost": 0 }

    if budget == 0:
        return output

    tmp_file = "input/tmp-G.csv"
    # partial from left
    G_L = np.copy(rewards)
    np.savetxt(tmp_file, G_L, delimiter=",", fmt='%.0f')
    output_sc_L = opt_partial_row_single_column_novertical_cpp(m, n, budget)

    G_R = np.copy(rewards)
    G_R = np.flip(G_R, 1)
    np.savetxt(tmp_file, G_R, delimiter=",", fmt='%.0f')
    output_sc_R = opt_partial_row_single_column_novertical_cpp(m, n, budget)

    js_L = output_sc_L["js"]
    js_R = output_sc_R["js"]

    full_rows = []
    full_rows = np.array(full_rows, dtype=int)

    # check at least two "full-rows" (always due to the while...)
    denom = 2
    len_full_rows = len(full_rows)
    while len_full_rows == 0:
        max_permitted_row = min(m, int((budget-2*n+4)/2))
        for i in range(0, max_permitted_row):
            if js_L[i] + js_R[i] > ((n+2)/denom):
                full_rows = np.append(full_rows, i)

        val_selected_rows = np.zeros(len(full_rows))
        for i in range(0, len(full_rows)):
            sel_i = full_rows[i]
            val_selected_rows[i] = np.sum(rewards[sel_i])

        if len(full_rows) % 2 == 1:
            worst_row = np.argsort(val_selected_rows)[0]

            full_rows = np.delete(full_rows, worst_row)
            val_selected_rows = np.delete(val_selected_rows, worst_row)

        len_full_rows = len(full_rows)
        denom = denom+1

    while True:
        traversed = len(full_rows)
        last_row = max(full_rows)
        reward_1 = np.sum(val_selected_rows)
        cost_1 = 2 * last_row + traversed * (n - 1)

        if cost_1 > budget:
            full_rows = full_rows[:-2]
            val_selected_rows = val_selected_rows[:-2]
        else:
            break

    budget_res = budget - cost_1
    if budget_res % 2 == 1:
        budget_res = budget_res - 1

    if budget_res > 0:
        S_t = np.sort(full_rows)

        # G_L
        tmp = np.copy(rewards[:, 0:(n / 2)])
        G_L = np.zeros((last_row - traversed + 1, np.size(tmp, 1)))
        tmp_i = 0
        for i in range(0, last_row):
            if i not in S_t:
                G_L[tmp_i] = tmp[i]
                tmp_i = tmp_i + 1

        # G_L "lower left corner"
        needed_B = (m - last_row - 1) * 2
        cost_2 = budget_res
        if budget_res >= needed_B:
            # it is not free for "go and back" there...
            budget_res = budget_res - needed_B
            G_L_plus = np.zeros((m - last_row - 1, np.size(tmp, 1)))
            tmp_i = 0
            for i in range(last_row + 1, m):
                G_L_plus[tmp_i] = tmp[i]
                tmp_i = tmp_i + 1

            G_L = np.block([[G_L], [G_L_plus]])

        # G_R
        tmp = np.copy(rewards[:, (n / 2):n])
        G_R = []
        for i in range(0, traversed - 1, 2):
            for k in range(S_t[i] + 1, S_t[i + 1]):
                G_R.append(tmp[k])
        G_R = np.array(G_R)

        if len(G_R > 0):
            G_R = np.flip(G_R, 1)
            G = np.block([[G_L], [G_R]])
        else:
            G = np.copy(G_L)

        if len(G > 0):
            [m_G, n_G] = G.shape

            tmp_file = "input/tmp-G.csv"
            np.savetxt(tmp_file, G, delimiter=",", fmt='%.0f')

            fast = True
            if fast:
                output_sc = opt_partial_row_single_column_novertical_cpp(m_G, n_G, budget_res)
            else:
                output_sc = opt_partial_row_single_column_novertical(m_G, n_G, budget_res)

            reward_2 = output_sc["reward"]
            # cost_2 = output_sc["cost"]

            cost_tot = cost_1 + cost_2
            reward_tot = reward_1 + reward_2

            output["reward"] = reward_tot
            output["cost"] = cost_tot

    return output


# {osc_L + osc_R}\{1,m} => osc
def heuristic_3(input, budget):
    rewards = read_input_from_file_id(input)
    [m, n] = rewards.shape
    output = {"cycle": [], "reward": 0, "cost": 0 }

    max_permitted_row = min(m, int((budget - 2 * n) / 2))

    if max_permitted_row < 0:
        return output

    reward_1 = np.sum(rewards[0]) + np.sum(rewards[max_permitted_row-1])
    cost_1 = 2*(n-1) + 2*(max_permitted_row-1)

    tmp = rewards[1:max_permitted_row-1]
    # G_L + G_R
    G1 = tmp[:, :(n/2)]
    G2 = np.flip(np.flip(tmp[:, -(n/2):], 1), 0)
    G = np.block([[G1], [G2]])

    budget_res = budget - cost_1

    if len(G > 0):
        [m_G, n_G] = G.shape

        tmp_file = "input/tmp-G.csv"
        np.savetxt(tmp_file, G, delimiter=",", fmt='%.0f')

        fast = True
        if fast:
            output_sc = opt_partial_row_single_column_novertical_cpp(m_G, n_G, budget_res)
        else:
            output_sc = opt_partial_row_single_column_novertical(m_G, n_G, budget_res)

        reward_2 = output_sc["reward"]
        # cost_2 = output_sc["cost"]

        cost_tot = cost_1 + budget_res
        reward_tot = reward_1 + reward_2

        output["reward"] = reward_tot
        output["cost"] = cost_tot

    return output


# ofr
def heuristic_4(input, budget):
    return opt_full_row(input, budget)


# osc
def heuristic_5(input, budget):
    return opt_partial_row_single_column_cpp(input, budget)


def opt_partial_row_single_column_novertical_cpp(m, n, budget):
    cmd = "oprsc.exe 2 0 %d %d %d" % (budget, m, n)
    output_sb = subprocess.check_output(cmd, shell=True)
    tmp_file = "input/tmp-G.csv"
    os.remove(tmp_file)
    cost = float(output_sb.strip().split(", ")[1].split("=")[1])
    reward = float(output_sb.strip().split(", ")[0].split("=")[1])

    js_str = output_sb.strip().split(", ")[2].split("=")[1]
    js_str = js_str[0:len(js_str)-1]
    js = np.array(js_str.split(","), dtype=int)

    output_oprsc = {"js": js, "reward": reward, "cost": cost}
    return output_oprsc


def opt_partial_row_single_column_novertical(m, n, budget):
    reader = csv.reader(open("input/tmp-G.csv", "rb"), delimiter=",")
    x = list(reader)
    rewards = np.array(x).astype("float")
    r = rewards
    B = budget

    output = {
        "cycle": [],
        "js": [],
        "reward": 0,
        "cost": 0,
    }

    V = np.empty((0, 2), int)
    V = np.append(V, [[0, 0]], axis=0)
    js = []

    T = np.cumsum(rewards, 1)

    nR = B/2+1
    R = np.zeros((m, nR), dtype=int)
    S = np.zeros((m, nR), dtype=int)

    # first row of R
    for b in range(0, nR):
        if b < n:
            R[0][b] = T[0][b]
            S[0][b] = b
        else:
            R[0][b] = T[0][n-1]
            S[0][b] = n-1

    for i in range(1, m):
        for b in range(0, nR):
            max_v = -100
            arg_v = -1
            for j in range(0, n):
                idx = b - j
                if idx >= 0:
                    if R[i - 1][idx] != -1:
                        v = R[i - 1][idx] + T[i][j]
                        if v > max_v:
                            max_v = v
                            arg_v = j

            R[i][b] = max_v
            S[i][b] = arg_v

    # print(R)
    # print(S)

    reward = R[m-1][nR-1]
    cost = 0
    j = nR-1
    # Sj = np.zeros(m)
    for i in range(m-1, -1, -1):
        m_j = S[i][j]
        # Sj[i] = m_j
        cost = cost + 2*m_j
        js.append(m_j)

        if m_j != 0:
            V = np.append(V, [[i, 0]], axis=0)
            V = np.append(V, [[i, m_j]], axis=0)
            V = np.append(V, [[i, 0]], axis=0)
        else:
            V = np.append(V, [[i, 0]], axis=0)

        j = j - S[i][j]

    # print("c=%d" % cost)
    output["cycle"] = V
    output["js"] = js
    output["reward"] = reward
    output["cost"] = cost
    return output


def opt_full_row(input, budget):
    rewards = read_input_from_file_id(input)
    [m, n] = rewards.shape
    output = {"cycle": [], "traversed": [], "reward": 0, "cost": 0}

    V = np.empty((0, 2), int)
    # T = 0
    v = [0, 0]

    R = np.zeros(m)
    for i in range(0, m):
        R[i] = np.sum(rewards[i])

    # rows belonging to the solution 0 no, 1 yes [i-th solution]
    S = np.zeros((m, m))

    for mp in range(0, m):
        Bp = budget - 2*mp
        k = int(floor(Bp / (n-1)))
        if k % 2 == 1:
            k = k - 1

        if k <= 0:
            continue

        if k > mp:
            for i in range(0, mp+1):
                S[mp][i] = 1
            continue

        Rt = R[0:mp]
        # Rt_s = np.sort(Rt)[::-1]
        Rt_sarg = np.argsort(Rt)[::-1]
        for i in range(0, k-1):
            S[mp][Rt_sarg[i]] = 1
        S[mp][mp] = 1

    # print(S)
    collected = np.zeros(m)
    for s in range(0, m):
        for i in range(0, m):
            if S[s][i] == 1:
                collected[s] = collected[s] + np.sum(rewards[i])

    # print(collected)

    best_sol = max(collected)
    best_arg = np.argmax(collected)

    if best_arg == 0:
        output["cycle"] = V
        output["traversed"] = []
        output["reward"] = 0
        output["cost"] = 0
        return output

    rows_sol = int(np.sum(S[best_arg]))
    if rows_sol % 2 == 1:
        rows_sol = rows_sol + 1
    cost = (n - 1)*rows_sol + 2*best_arg

    reward = best_sol

    V = np.append(V, [v], axis=0)
    for i in range(0, m):
        if S[best_arg][i] == 1:
            v[0] = i
            V = np.append(V, [v], axis=0)
            if v[1] == 0:
                v[1] = n - 1
            else:
                v[1] = 0
            V = np.append(V, [v], axis=0)
    V = np.append(V, [[0, 0]], axis=0)

    # clean
    ok = True
    while ok:
        i = 0
        while i < len(V) - 1:
            if (V[i][0] == V[i + 1][0]) and (V[i][1] == V[i + 1][1]):
                V = np.delete(V, i, 0)
                break
            i = i + 1
        if i == len(V) - 1:
            ok = False

    sol_size = int(np.sum(S[best_arg]))
    Rt = R[0:best_arg+1]
    # Rt_s = np.sort(Rt)[::-1]
    Rt_sarg = np.argsort(Rt)[::-1]
    s = Rt_sarg[0:sol_size]

    output["cycle"] = V
    output["traversed"] = s
    output["reward"] = reward
    output["cost"] = cost
    return output


def opt_partial_row_single_column_cpp(input, budget):
    rewards = read_input_from_file_id(input)
    [m, n] = rewards.shape
    output = {"cycle": [], "traversed": [], "reward": 0, "cost": 0}

    cmd = "oprsc.exe 1 %d %d %d %d" % (input, budget, m, n)
    output_sb = subprocess.check_output(cmd, shell=True)
    cost = float(output_sb.strip().split(", ")[1].split("=")[1])
    reward = float(output_sb.strip().split(", ")[0].split("=")[1])
    output_oprsc = {
        "reward": reward,
        "cost": cost
    }

    return output_oprsc


def opt_partial_row_single_column(input, budget):
    rewards = read_input_from_file_id(input)
    [m, n] = rewards.shape
    output = {"cycle": [], "traversed": [], "reward": 0, "cost": 0}

    V = np.empty((0, 2), int)
    V = np.append(V, [[0, 0]], axis=0)

    T = np.cumsum(rewards, 1)

    nR = B/2+1
    R = np.zeros((m, nR), dtype=int)
    S = np.zeros((m, nR), dtype=int)

    # first row of R
    for b in range(0, nR):
        if b < n:
            R[0][b] = T[0][b]
            S[0][b] = b
        else:
            R[0][b] = T[0][n-1]
            S[0][b] = n-1

    for i in range(1, m):
        # print(i)
        for b in range(0, nR):
            if b < i:
                R[i][b] = -1
                S[i][b] = -1
            elif b == i:
                R[i][b] = 0
                S[i][b] = 0
            else:
                max_v = -100
                arg_v = -1
                for j in range(0, n):
                    idx = b-j-1
                    if idx >= 0:
                        if R[i-1][idx] != -1:
                            v = R[i-1][idx] + T[i][j]
                            if v > max_v:
                                max_v = v
                                arg_v = j

                R[i][b] = max_v
                S[i][b] = arg_v

    # print(R)
    # print(S)

    reward = np.max(R[:, np.size(R, axis=1)-1])
    max_m = np.argmax(R[:, np.size(R, axis=1)-1])

    cost = 2*max_m
    j = nR-1
    for i in range(max_m, -1, -1):
        m_j = S[i][j]
        cost = cost + 2*m_j

        if m_j != 0:
            V = np.append(V, [[i, 0]], axis=0)
            V = np.append(V, [[i, m_j]], axis=0)
            V = np.append(V, [[i, 0]], axis=0)
        else:
            V = np.append(V, [[i, 0]], axis=0)

        j = j - (S[i][j]+1)

    # print("c=%d" % cost)
    output["cycle"] = V
    output["reward"] = reward
    output["cost"] = cost
    return output


def is_vertex_feasible(i, j, T, B, v, m, n):
    cost = 0

    # vertical from v to the i-th row
    cost = cost + fabs(i - v[0])

    if v[1] == 0:
        # if v is in the left side
        # traverse the row up to the j-th column (and come back)
        cost = cost + 2*j
    else:
        # if v is in the right side
        # traverse the entire row
        cost = cost + n - 1

    # vertical from current to home
    cost = cost + i

    # print("v=[%d,%d], (%d,%d) cost=%d, T=%d, B=%d" % (v[0], v[1], i, j, cost, T, B))

    if T + cost > B:
        return 0
    else:
        return 1


def greedy_partial_row_single_column(input, budget):
    rewards = read_input_from_file_id(input)
    [m, n] = rewards.shape
    output = {"cycle": [], "traversed": [], "reward": 0, "cost": 0}

    V = np.empty((0, 2), int)
    T = 0
    v = [0, 0]
    V = np.append(V, [v], axis=0)

    L = np.zeros((m, n))
    for i in range(0, m):
        L[i][0] = r[i][0]
        for j in range(1, n):
            L[i][j] = L[i][j-1] + r[i][j]

    Fij = np.ones((m, n))
    # Fr = np.ones(m)

    Sij = np.zeros((m, n))
    # Sr = np.zeros(m)

    while np.sum(Fij) > 0:
        for i in range(0, m):
            for j in range(0, n):
                if Sij[i][j] == 0:
                    Fij[i][j] = is_vertex_feasible(i, j, T, B, v, m, n)

        if np.sum(Fij) == 0:
            break

        Lp = np.ones((m, n)) * -1
        for i in range(0, m):
            for j in range(0, n):
                if Fij[i][j] == 1 and Sij[i][j] == 0:
                    c = fabs(i - v[0]) + 2*j
                    if c != 0:
                        Lp[i][j] = L[i][j] / c

        next = np.unravel_index(np.nanargmax(Lp), np.array(Lp).shape)

        # from v to the next-th row
        T = T + fabs(next[0] - v[0])
        v[0] = next[0]
        V = np.append(V, [v], axis=0)

        v[1] = next[1]
        V = np.append(V, [v], axis=0)
        v[1] = 0
        V = np.append(V, [v], axis=0)

        T = T + 2 * fabs(next[1] - v[1])

        for j in range(0, next[1] + 1):
            Fij[next[0]][j] = 0
            Sij[next[0]][j] = 1

        # update the cumulative rewards
        mask = ~np.array(Sij, dtype=bool)
        tmp = r.copy()
        tmp[~mask] = 0

        for i in range(0, m):
            L[i][0] = tmp[i][0]
            for j in range(1, n):
                L[i][j] = L[i][j - 1] + tmp[i][j]

    if (v[0] != 0) or (v[1] != 0):
        V = np.append(V, [[0, 0]], axis=0)
        T = T + next[0]

    # clean
    ok = True
    while ok:
        i = 0
        while i < len(V)-1:
            if (V[i][0] == V[i+1][0]) and (V[i][1] == V[i+1][1]):
                V = np.delete(V, i, 0)
                break
            i = i+1
        if i == len(V)-1:
            ok = False

    # fix
    max_i = -1
    for i in range(0, m):
        if Sij[i][0] == 1:
            max_i = i
    for i in range(0, max_i):
        Sij[i][0] = 1

    cost = T
    reward = 0
    for i in range(0, m):
        for j in range(0, n):
            if Sij[i][j] == 1:
                reward = reward + r[i][j]

    output["cycle"] = V
    output["reward"] = reward
    output["cost"] = cost
    return output


def execute_batch(function_name, input):
    rewards = read_input_from_file_id(input)
    [m, n] = rewards.shape

    min_budget = 2 * n
    step_budget = n
    max_budget = n * m + 3 * m

    # budgets = range(min_budget, max_budget+step_budget, step_budget)
    budgets = np.linspace(0, max_budget, num=20)
    size = len(budgets)

    rewards = np.zeros(size)
    costs = np.zeros(size)

    if function_name == 'gfr' or function_name == 'gpr':
        _call_matlab_function(function_name, input, min_budget, step_budget, max_budget)
        return

    for b in range(0, size):
        print("%d/%d" % (b+1, size))

        budget = int(budgets[b])

        if function_name == 'ofr':
            output = opt_full_row(input, budget)
        elif function_name == 'oprsc':
            output = opt_partial_row_single_column_cpp(input, budget)
        elif function_name == 'gprsc':
            output = greedy_partial_row_single_column(input, budget)
        elif function_name == 'hprgc':
            output = heuristic_partial_row(input, budget)

        rewards[b] = output["reward"]
        costs[b] = output["cost"]

    file_path = 'output/%s/out_%s-%d.txt' % (function_name, function_name, input)
    f = open(file_path, "w")
    for i in range(0, len(rewards)):
        f.write("%.2f,%.2f,%.2f\n" % (int(budgets[i]), costs[i], rewards[i]))
    f.close()


def greedy_full_row(input, budget):
    h = win32com.client.Dispatch('matlab.application')
    # h.Execute("addpath('C:\Users\\franc\Desktop\droneranging\papers\ICRA 2020\\robots');")

    func = "greedy_full_row(input, budget)"
    func = func.replace("input", str(input))
    func = func.replace("budget", str(budget))

    out = h.Execute(func)
    string = str(out).strip().split("|")[1].split(",")
    output = {
        "cycle": [],
        "reward": float(string[2].split("=")[1]),
        "cost": int(string[1].split("=")[1]),
    }
    return output


def greedy_partial_row(input, budget):
    h = win32com.client.Dispatch('matlab.application')
    # h.Execute("addpath('C:\Users\\franc\Desktop\droneranging\papers\ICRA 2020\\robots');")

    func = "greedy_partial_row(input, budget)"
    func = func.replace("input", str(input))
    func = func.replace("budget", str(budget))

    out = h.Execute(func)
    string = str(out).strip().split("|")[1].split(",")
    output = {
        "cycle": [],
        "reward": float(string[2].split("=")[1]),
        "cost": int(string[1].split("=")[1]),
    }
    return output


def _call_matlab_function(name, input, min_budget, step_budget, max_budget):
    h = win32com.client.Dispatch('matlab.application')
    # h.Execute("addpath('C:\Users\\franc\Desktop\droneranging\papers\ICRA 2020\\robots');")

    if name == 'gfr':
        function_name = 'greedy_full_row_batch'
    elif name == 'gpr':
        function_name = 'greedy_partial_row_batch'

    func = "function_name(input, min_budget, step_budget, max_budget)"
    func = func.replace("function_name", str(function_name))
    func = func.replace("input", str(input))
    func = func.replace("min_budget", str(min_budget))
    func = func.replace("step_budget", str(step_budget))
    func = func.replace("max_budget", str(max_budget))

    h.Execute(func)

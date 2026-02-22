
input = 20;
budget = 190;


[budgets, costs, rewards] = greedy_full_row(input, budget);
fprintf('FULL -> budget=%d, cost=%d, reward=%d\n', budgets, costs, rewards);

[budgets, costs, rewards] = greedy_partial_row(input, budget);
fprintf('PARTIAL -> budget=%d, cost=%d, reward=%d\n', budgets, costs, rewards);

greedy_full_row_batch(input, 0, 0, 12000)
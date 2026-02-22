function greedy_partial_row_batch(input, min_budget, step_budget, max_budget)

    rewards_out = [];
    costs_out = [];
    budgets_out = [];
    
    budgets = linspace(0, max_budget, 20);
    
%     for budget=min_budget:step_budget:max_budget+step_budget
    for i=1:size(budgets,2)
        budget = floor(budgets(i));
        
        [b, c, r] = greedy_partial_row(input, budget);
        
        rewards_out = [rewards_out, r];
        costs_out = [costs_out, c];
        budgets_out = [budgets_out, b];
    end
    
    
    out_filename = sprintf('C:\\Users\\franc\\Desktop\\droneranging\\papers\\ICRA 2020\\robots\\output\\gpr\\out_gpr-%d.txt', input);
    fileID = fopen(out_filename,'w');
    for i=1:size(rewards_out, 2)
        fprintf(fileID, '%.2f,%.2f,%.2f\n', budgets_out(i), costs_out(i), rewards_out(i));
    end
    
    fclose(fileID);
end

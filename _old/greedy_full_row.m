function [budget_out, cost_out, reward_out] = greedy_full_row(input, budget)

    reward_out = 0;
    cost_out = 0;
    budget_out = 0;
    
    input_filename = sprintf('C:\\Users\\franc\\Desktop\\droneranging\\papers\\ICRA 2020\\robots\\input\\input-%d.csv', input); 
    reward_map = csvread(input_filename);

    %% Compile reward_out and cost_out for each row
    beginning_row = ceil(1/size(reward_map,2));
    ending_row = ceil(1/size(reward_map,2));
    num_vines_per_row = size(reward_map, 2);
    size_row = size(reward_map, 1); 
    size_column = size(reward_map, 2);
    for i=1:size(reward_map,1)
        compiled_reward_out(i) = sum(reward_map(i,:));
    end
    compiled_row_cost = num_vines_per_row * 1;
    [reward_out, index] = sort(compiled_reward_out,'descend');
    reward_out = [reward_out', index'];
    current_row = beginning_row;
    current_side = 1;
    total_cost = 0;
    total_reward = 0;
    tour = [];
    while length(reward_out) > 0
        if current_side == 1
            loop_cost = 1*abs(current_row - reward_out(1, 2)) + 2*compiled_row_cost + 1*abs(reward_out(1,2) - ending_row);
            if loop_cost <= (budget - total_cost)
                total_cost = total_cost + 1*abs(current_row - reward_out(1, 2)) + compiled_row_cost;
                total_reward = total_reward + reward_out(1, 1);
                tour = [tour; reward_out(1, 2)];
                current_row = reward_out(1, 2);
                current_side = 2;
                reward_out(1, :) = [];
            else
                reward_out(1, :) = [];
            end
        elseif current_side == 2
            loop_cost = 1*abs(current_row - reward_out(1, 2)) + compiled_row_cost + 1*abs(reward_out(1,2) - ending_row);
            if loop_cost <= (budget - total_cost)
                total_cost = total_cost + 1*abs(current_row - reward_out(1, 2)) + compiled_row_cost;
                total_reward = total_reward + reward_out(1, 1);
                tour = [tour; reward_out(1, 2)];
                current_row = reward_out(1, 2);
                current_side = 1;
                reward_out(1, :) = [];
            else
                reward_out(1, :) = [];
            end
        end
    end
    if current_side == 2 % Will only be at side 2 if the row needed to reach the 1 has already been traversed
        total_cost = total_cost + 1*abs(current_row - ending_row) + compiled_row_cost;
        %tour = [tour; 1];
    end
    if current_side == 1
        total_cost = total_cost + 1*abs(current_row - ending_row);
    end

    now = 1;
    reward = 0;
    real_tour = [];
    cost = 0;
    for i=1:length(tour)
        now_row = ceil(now/size(reward_map,2));
        now_col = rem((now-1),size(reward_map,2)) + 1;
        if now_row < tour(i)
            for j=now_row:(tour(i)-1)
                real_tour = [real_tour; (j-1)*size(reward_map,2) + now_col];
                reward = reward + reward_map(j, now_col);
                reward_map(j, now_col) = 0;
                cost = cost + 1;
            end
        else
            for j=now_row:-1:(tour(i)+1)
                real_tour = [real_tour; (j-1)*size(reward_map,2) + now_col];
                reward = reward + reward_map(j, now_col);
                reward_map(j, now_col) = 0;
                cost = cost + 1;
            end
        end
        if now_col == 1
            for j=now_col:(size(reward_map,2)-1)
                real_tour = [real_tour; (tour(i)-1)*size(reward_map,2) + j];
                reward = reward + reward_map(tour(i), j);
                reward_map(tour(i), j) = 0;
                cost = cost + 1;
            end
            now = (tour(i)-1)*size(reward_map,2) + size(reward_map,2);
        else
            for j=now_col:-1:2
                real_tour = [real_tour; (tour(i)-1)*size(reward_map,2) + j];
                reward = reward + reward_map(tour(i), j);
                reward_map(tour(i), j) = 0;
                cost = cost + 1;
            end
            now = (tour(i)-1)*size(reward_map,2) + 1;
        end
    end
    now_row = ceil(now/size(reward_map,2));
    now_col = rem((now-1),size(reward_map,2)) + 1;
    if now_row < ending_row
        for j=now_row:(ending_row)
            real_tour = [real_tour; (j-1)*size(reward_map,2) + now_col];
            reward = reward + reward_map(j, now_col);
            reward_map(j, now_col) = 0;
            cost = cost + 1;
        end
    else
        for j=now_row:-1:(ending_row)
            real_tour = [real_tour; (j-1)*size(reward_map,2) + now_col];
            reward = reward + reward_map(j, now_col);
            reward_map(j, now_col) = 0;
            cost = cost + 1;
        end
    end
    now = real_tour(end);
    now_row = ceil(now/size(reward_map,2));
    now_col = rem((now-1),size(reward_map,2)) + 1;
    if now_col ~= 1
        for j=now_col-1:-1:1
            real_tour = [real_tour; (now_row-1)*size(reward_map,2) + j];
            reward = reward + reward_map(now_row, j);
            reward_map(now_row, j) = 0;
            cost = cost + 1;
        end
    end

    %% Eliminate duplicates
    i = 2;
    while i <= length(tour)
        if tour(i-1) == tour(i)
            tour(i-1) = [];
        else
            i = i + 1;
        end
    end
    
    
    reward_out =  total_reward;
    cost_out = total_cost;
    budget_out = budget;
    fprintf('|budget=%d, cost=%d, reward=%d|\n', budget_out, cost_out, reward_out);
end

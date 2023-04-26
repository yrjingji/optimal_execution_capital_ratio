function all_combinations = get_combinations(vectors)
    n = numel(vectors);
    all_combinations = [];

    function traverse(index, temp_combination)
        if index > n
            all_combinations(end + 1, :) = temp_combination;
        else
            for i = 1:numel(vectors{index})
                temp_combination(index) = vectors{index}(i);
                traverse(index + 1, temp_combination);
            end
        end
    end

    traverse(1, zeros(1, n));
end


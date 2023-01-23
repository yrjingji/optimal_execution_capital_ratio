for k=0:0.01:1
    real = check_naive_single_feasible(k);
    if any(real <0.95)
        j = k;
    end
end
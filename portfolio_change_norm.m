sigma = [0.9 0.4;0.4 0.9];
% we want to maintatin that the original variance of the portfolio is
% constant with different correlations given
X_0 = [305,300];
k = 0.2;
cons_variance = X_0*sigma*X_0';
rho_seq = [-0.99:0.01:-0.85,-0.84:0.1:1];
parfor j = 1:length(rho_seq)
    %fix one variance to be constant and change the other one.
    correaltion_matrix = [1 rho_seq(j); rho_seq(j) 1];
    var_same = cons_variance/(X_0*correaltion_matrix*X_0');
    sigma = var_same * correaltion_matrix;
    [opt_val_portfolio_rho(j),sol_first_rho(j,:),sol_second_rho(j,:)] = ...
        general_model_portfolio_5time(k*ones(1,2),sigma);
    [opt_special_portfolio_rho(j), sol_first_special_rho(j,:), sol_second_special_rho(j,:)]=...
        special_trading_portfolio_2asset_5time(k*ones(1,2),sigma);
    [naive_feasible(j)] = all(check_naive_portfolio(k*ones(1,2),sigma)>0.95);
    if naive_feasible(j) ==1
        naive_value(j) = naive(sigma);
    else
        naive_value(j) = -Inf;
    end
end
plot([-0.99:0.01:-0.85,-0.84:0.05:0.5],opt_val_portfolio_rho,'b');
hold on
plot([-0.99:0.01:-0.85,-0.84:0.05:0.5],opt_special_portfolio_rho,'m');
plot([-0.99:0.01:-0.85,-0.84:0.05:0.5],naive_value,'r');
legend('optimal trading strategy','affine trading strategy','naive trading strategy');
xlabel('correlation');
ylabel('expected revenue');
title('portfolio')
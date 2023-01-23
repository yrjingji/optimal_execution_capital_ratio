j=1;
for mu=0:0.01:0.5
    [opt_single_meanvar_val_noconstraint(j),opt_single_meanvar_mean_noconstraint(j), ...
        opt_single_meanvar_variance_noconstraint(j),opt_single_meanvar_strategy_noconstraint(j,:)] = ...
        direct_chance_5time_mean_var_noconstraint(mu);
    j = j+1;
end
initial_price = 17;
%standard deviation
sigma = 0.7;
tau = 1;
l = 1200;
X_0 = 305;
S_0 = 250;
beta = .09;
gamma = .072;
%objective function
negative_P0 = 0.5*gamma*ones(5,5) + diag(ones(1,5)*(beta - 0.5*gamma));
negative_q0 = -initial_price*ones(1,5);
naive_trading = ones(5,1)*S_0/5;
Q(1,:) = sqrt(tau)*sigma*[0 1 1 1 1];
Q(2,:) = sqrt(tau)*sigma*[0 0 1 1 1];
Q(3,:) = sqrt(tau)*sigma*[0 0 0 1 1];
Q(4,:) = sqrt(tau)*sigma*[0 0 0 0 1];
var_coeff = Q*naive_trading;
mean =  -quad_form(naive_trading,negative_P0) - dot(negative_q0,naive_trading);
variance = var_coeff'*var_coeff;
%opt_sharpe_ratio = opt_single_meanvar_mean_noconstraint./sqrt(opt_single_meanvar_variance_noconstraint);
%check feasibility of solutions:
for j=1:51
    feasible_not(j) = any(check_single_feasible(opt_single_meanvar_strategy_noconstraint(j,:)',0.5) >=0.95);
end
 plot(0:0.01:0.21,opt_single_meanvar_mean_noconstraint(1:22)./sqrt(opt_single_meanvar_variance_noconstraint(1:22)),'-b');
hold on
 plot(0:0.01:0.21,opt_sharpe_ratio(1:22),'-r');
 mu_seq= 0:0.01:0.21;
xline(mu_seq(find(feasible_not ==1,1,'last')),'--m','cut-off point');

% ylabel('objective value')
%legend('optimal trading strategy without chance constraint','optimal trading strategy with chance constraint');
% xlabel('risk aversion')
%title('single asset')
yline( mean/sqrt(variance),'-k');
ylabel('mean/std')
legend('optimal trading strategy without chance constraint','optimal trading strategy with chance constraint','threshold','naive trading strategy');
xlabel('risk aversion')
xlim([0,0.24])
title('single asset')
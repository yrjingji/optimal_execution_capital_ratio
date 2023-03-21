mu_seq= 0:0.01:0.21;
parfor j=1:length(mu_seq)
    [opt_single_meanvar_val(j),opt_single_meanvar_mean(j), ...
        opt_single_meanvar_variance(j),opt_single_meanvar_strategy(j,:),single_prob_split(j,:)] = ...
        direct_chance_5time_mean_var(0.3,mu_seq(j));
    [opt_single_meanvar_val_noconstraint(j),opt_single_meanvar_mean_noconstraint(j), ...
        opt_single_meanvar_variance_noconstraint(j),opt_single_meanvar_strategy_noconstraint(j,:)] = ...
        direct_chance_5time_mean_var_noconstraint(mu_seq(j));
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
for j=1:length(mu_seq)
    feasible_not(j) = any(check_single_feasible(opt_single_meanvar_strategy_noconstraint(j,:)',0.5) >=0.95);
end
opt_sharpe_ratio_cons =  opt_single_meanvar_mean./sqrt(opt_single_meanvar_variance);
opt_sharpe_ratio_noncons = opt_single_meanvar_mean_noconstraint./sqrt(opt_single_meanvar_variance_noconstraint);
naive_sharpe = mean/sqrt(variance);
%---- plot the sharpe ratio
figure(1)
plot(mu_seq,opt_sharpe_ratio_noncons,'-r');
hold on
plot(mu_seq,opt_sharpe_ratio_cons,'-b');
yline(naive_sharpe,'-k')
xline(mu_seq(find(feasible_not ==1,1,'last')),'--m','cut-off point');
xlim([0,0.21]);
ylabel('mean/std')
legend('optimal trading strategy without chance constraint','optimal trading strategy with chance constraint','naive trading strategy','threshold');
xlabel('risk aversion')
title('single asset')


%-------objective value generator
figure(2)
plot(0:0.01:0.21,opt_single_meanvar_val_noconstraint,'-r');
hold on
plot(mu_seq,opt_single_meanvar_val,'-b');
plot(0:0.01:0.21,mean-mu_seq*variance,'-k')
xline(mu_seq(find(feasible_not ==1,1,'last')),'--m','cut-off point');
xlim([0,0.21]);
ylabel('objective value')
legend('optimal trading strategy without chance constraint','optimal trading strategy with chance constraint','naive trading strategy','threshold');
xlabel('risk aversion')
title('single asset: max mean-var')

j=1;
for mu=0:0.01:1
    [opt_single_meanvar_val(j),opt_single_meanvar_mean(j), ...
        opt_single_meanvar_variance(j),opt_single_meanvar_strategy(j,:)] = ...
        direct_chance_5time_mean_var_fixedalpha(0.5,mu,alpha_1,alpha_2,alpha_3,alpha_4);
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
% opt_sharpe_ratio = opt_single_meanvar_mean./sqrt(opt_single_meanvar_variance);
% plot(0:0.01:1,opt_sharpe_ratio,'-b');
% yline( mean/sqrt(variance),'-r');
% ylabel('mean/std')
% legend('optimal trading strategy','naive trading strategy');
% xlabel('risk aversion')
% title('single asset')



%-------
mu_seq= 0:0.01:0.21;
plot(0:0.01:0.21,opt_single_meanvar_val(1:22),'-b');
hold on
plot(0:0.01:0.21,opt_single_meanvar_val_noconstraint(1:22),'-r');
plot(0:0.01:0.21,mean-mu_seq*variance,'-k')
xline(mu_seq(find(feasible_not ==1,1,'last')),'--m','cut-off point');
%ylim([-7000,2500])
xlim([0,0.21]);
ylabel('objective value')
legend('optimal trading strategy without chance constraint','optimal trading strategy with chance constraint','threshold','naive trading strategy');
xlabel('risk aversion')
title('single asset')

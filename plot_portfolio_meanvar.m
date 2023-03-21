j=1;
sigma = [0.9 0.2;0.2 0.4];
for mu=0:0.01:0.42
    [opt_portfolio_meanvar_val_noconstraint(j),opt_portfolio_meanvar_mean_noconstraint(j), ...
        opt_portfolio_meanvar_variance_noconstraint(j),opt_portfolio_meanvar_first_noconstraint(j,:),opt_portfolio_meanvar_second_noconstraint(j,:)] = ...
        general_model_portfolio_5time_meanvar_noconstraint(sigma,mu);
    j = j+1;
end
% plot(0:0.01:0.42,opt_portfolio_meanvar_val_noconstraint(1:43),'-b');
% hold on
% plot(0:0.01:0.42,opt_portfolio_meanvar_val,'-k'); 
% plot(0:0.01:0.42,opt_portfolio_meanvar_val_special,'--r')
% mu_seq= 0:0.01:0.42;
% xline(mu_seq(find(feasible_not ==1,1,'last')),'--m','cut-off point');
% ylabel('objective value');
% legend('optimal trading strategy without constraint','optimal trading strategy with constraint','affine trading strategy');
% xlabel('risk aversion')
% title('portfolio')

plot(0:0.01:0.42,(opt_portfolio_meanvar_val-opt_portfolio_meanvar_val_noconstraint(1:43))./opt_portfolio_meanvar_val_noconstraint(1:43),'-b');
hold on
plot(0:0.01:0.42,(opt_portfolio_meanvar_val_special-opt_portfolio_meanvar_val_noconstraint(1:43))./opt_portfolio_meanvar_val_noconstraint(1:43),'-r')
mu_seq= 0:0.01:0.42;
xline(mu_seq(find(feasible_not ==1,1,'last')),'--k','cut-off point');
ylabel('objective value');
legend('optimal trading strategy without constraint','optimal trading strategy with constraint','affine trading strategy');
xlabel('risk aversion')
title('portfolio')



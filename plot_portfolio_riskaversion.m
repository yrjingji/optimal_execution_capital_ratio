sigma = [0.9 0.2;0.2 0.9];
%capital ratio requirement changed here
k = 0.3;
mu_seq = [0:0.01:0.43];
parfor j=1:length(mu_seq)
    [opt_portfolio_meanvar_val(j),opt_portfolio_meanvar_mean(j), ...
        opt_portfolio_meanvar_variance(j),opt_portfolio_meanvar_first(j,:),opt_portfolio_meanvar_second(j,:),...
        prob_split_product_meanvar(j,:)] = ...
        general_model_portfolio_5time_meanvar(k*ones(1,2),sigma,mu_seq(j));
    [opt_portfolio_meanvar_val_special(j),opt_portfolio_meanvar_mean_special(j), ...
        opt_portfolio_meanvar_variance_special(j),opt_portfolio_meanvar_first_special(j,:),opt_portfolio_meanvar_second_special(j,:),...
        prob_split_product_meanvar_affine(j,:)]=...
        special_model_portfolio_5time_meanvar(k*ones(1,2),sigma,mu_seq(j));
    [opt_val_meanvar_noconstraint(j),~,~,sol_first_noconstraint(j,:),sol_second_noconstraint(j,:)]= general_model_portfolio_5time_meanvar_noconstraint(sigma,mu_seq(j));
end
initial_price = [17;17];
tau = 1;
l = 1200;
X_0 = [305;305];
S_0 = [250;200];
c1 = .07;
c2 = 0.04;
beta = c1*sigma;
gamma = c2*sigma;
%objective function (need modification)
for j=1:length(mu_seq)
    feasible_not(j) = any(check_portfolio_feasibility(sol_first_noconstraint(j,:),sol_second_noconstraint(j,:),k,sigma) >=0.95);
end
%-----plot objective value
figure(1)
plot(mu_seq,opt_val_meanvar_noconstraint,'-r');
hold on
plot(mu_seq,opt_portfolio_meanvar_val,'-b');
plot(mu_seq,opt_portfolio_meanvar_val_special,'-k')
xline(mu_seq(find(feasible_not ==1,1,'last')),'--m','cut-off point');
xlim([mu_seq(1),mu_seq(length(mu_seq))]);
ylabel('objective value')
legend('optimal trading strategy without chance constraint','optimal trading strategy with chance constraint','threshold');
xlabel('risk aversion')
title('portfolio')

%-----show the differenve among the value
figure(2)

plot(mu_seq,opt_portfolio_meanvar_val-opt_val_meanvar_noconstraint,'-b');
hold on
plot(mu_seq,opt_portfolio_meanvar_val_special-opt_val_meanvar_noconstraint,'-k')
xline(mu_seq(find(feasible_not ==1,1,'last')),'--m','cut-off point');
xlim([mu_seq(1),mu_seq(length(mu_seq))]);
ylabel('objective value')
legend('optimal trading strategy without chance constraint','optimal trading strategy with chance constraint','threshold');
xlabel('risk aversion')
title('portfolio')



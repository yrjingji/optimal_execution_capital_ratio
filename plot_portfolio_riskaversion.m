j=1;
sigma = [0.9 0.2;0.2 0.9];
U1 = chol(sigma);
k = 0.3;
mu_seq = [0:0.01:0.43];
parfor j=1:length(mu_seq)
%     [opt_portfolio_meanvar_val(j),opt_portfolio_meanvar_mean(j), ...
%         opt_portfolio_meanvar_variance(j),opt_portfolio_meanvar_first(j,:),opt_portfolio_meanvar_second(j,:)] = ...
%         general_model_portfolio_5time_meanvar(k*ones(1,2),sigma,mu);
    [opt_portfolio_meanvar_val_special(j),opt_portfolio_meanvar_mean_special(j), ...
        opt_portfolio_meanvar_variance_special(j),opt_portfolio_meanvar_first_special(j,:),opt_portfolio_meanvar_second_special(j,:)]=...
        special_model_portfolio_5time_meanvar(k*ones(1,2),sigma,mu_seq(j));
end
% initial_price = [17;17];
% tau = 1;
% l = 1200;
% X_0 = [305;305];
% S_0 = [250;200];
% c1 = .07;
% c2 = 0.04;
% beta = c1*sigma;
% gamma = c2*sigma;
% %objective function (need modification)
% A_obj = c1*sigma;
% B_obj = 0.5*c2*sigma;
% negative_P0 =[A_obj B_obj B_obj B_obj B_obj; B_obj A_obj B_obj B_obj B_obj;...
%     B_obj B_obj A_obj B_obj B_obj;B_obj B_obj B_obj A_obj B_obj;B_obj B_obj B_obj B_obj A_obj];
% negative_q0 = -[initial_price(1);initial_price(2);initial_price(1);initial_price(2);initial_price(1);initial_price(2);...
%     initial_price(1);initial_price(2);initial_price(1);initial_price(2)];
% s1 = S_0(1)/5*ones(1,5);
% s2 = S_0(2)/5*ones(1,5);
% trans_s = [s1(1),s2(1),s1(2),s2(2),s1(3),s2(3),s1(4),s2(4),s1(5),s2(5)];
% q1 = sqrt(tau)*U1'*(S_0-[s1(1);s2(1)]);
% q2 = sqrt(tau)*U1'*(S_0-[s1(1);s2(1)] -[s1(2);s2(2)]);
% q3 = sqrt(tau)*U1'*([s1(4);s2(4)]+[s1(5);s2(5)]);
% q4 = sqrt(tau)*U1'*([s1(5);s2(5)]);
% mean = -quad_form(trans_s,negative_P0) - dot(negative_q0,trans_s);
% variance = dot(q1,q1)+dot(q2,q2)+dot(q3,q3)+dot(q4,q4);
% opt_sharpe_ratio = opt_portfolio_meanvar_mean./sqrt(opt_portfolio_meanvar_variance);
% opt_sharpe_ratio_special = opt_portfolio_meanvar_mean_special./sqrt(opt_portfolio_meanvar_variance_special);
% plot(0:0.01:0.3,opt_sharpe_ratio,'-b');
% hold on
% plot(0:0.01:0.3,opt_sharpe_ratio_special,'-r');
% ylabel('mean/std')
% legend('optimal trading strategy','affine trading strategy');
% xlabel('risk aversion')
% 
% title('portfolio')
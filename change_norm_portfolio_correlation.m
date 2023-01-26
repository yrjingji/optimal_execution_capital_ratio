j=1;
sigma = [0.9 0.2;0.2 0.9];
% we want to maintatin that the original variance of the portfolio is
% constant with different correlations given
X_0 = [305,305];
k = 0.3;
cons_variance = X_0*sigma*X_0';
rho_seq =  -0.5:0.1:0.5;
parfor j = 1:length(rho_seq)
    cov = sqrt(sigma(1,1)*sigma(2,2))*rho_seq(j);
    sigma(1,2) = cov;
    sigma(2,1) = cov;
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
initial_price = [17,17];
tau = 1;
l = 1200;
S_0 = [250,200];
c1 = .07;
c2 = 0.04;
beta = c1*sigma;
gamma = c2*sigma;
%objective function
A_obj = c1*sigma;
B_obj = 0.5*c2*sigma;
negative_P0 =[A_obj B_obj B_obj B_obj B_obj; B_obj A_obj B_obj B_obj B_obj;...
    B_obj B_obj A_obj B_obj B_obj;B_obj B_obj B_obj A_obj B_obj;B_obj B_obj B_obj B_obj A_obj];
negative_q0 = -[initial_price(1);initial_price(2);initial_price(1);initial_price(2);initial_price(1);initial_price(2);...
    initial_price(1);initial_price(2);initial_price(1);initial_price(2)];
s1 = ones(5,1)*S_0(1)/5;
s2 = ones(5,1)*S_0(2)/5;
naive_trading = [s1(1),s2(1),s1(2),s2(2),s1(3),s2(3),s1(4),s2(4),s1(5),s2(5)];
value = -quad_form(naive_trading,negative_P0) - dot(negative_q0,naive_trading);
plot(-0.5:0.1:0.5,opt_val_portfolio_rho,'b');
hold on
plot(-0.5:0.1:0.5,opt_special_portfolio_rho,'m');
%yline(value,'r');
legend('optimal trading strategy','affine trading strategy');%,'naive trading strategy');
xlabel('correlation');
ylabel('expected revenue');
title('portfolio')
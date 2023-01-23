j=1;
sigma = [0.9,0.2;0.2 0.9];
for k=[0:0.05:0.65]
    [opt_val_portfolio(j),sol_first(j,:),sol_second(j,:)] = ...
        general_model_portfolio_5time(k*ones(1,2),sigma);
    [opt_special_portfolio(j), sol_first_special(j,:), sol_second_special(j,:)]=...
        special_trading_portfolio_2asset_5time(k*ones(1,2),sigma);
    j = j+1;
end
initial_price = [17,17];
tau = 1;
l = 1200;
X_0 = [305,305];
S_0 = [250,250];
c1 = .1;
c2 = 0.09;
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
plot([1.-[0:0.05:0.4,0.4:0.01:0.65]]/8,opt_val_portfolio,'b');
hold on
plot([1.-[0:0.05:0.4,0.4:0.01:0.65]]/8,opt_special_portfolio,'m');
legend('optimal trading strategy','affine trading strategy');
xlabel('requirement of capital ratio');
ylabel('expected revenue');
title('portfolio')
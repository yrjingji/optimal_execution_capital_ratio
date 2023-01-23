j=1;
for k=[0:0.05:0.4,0.4:0.01:0.68]
    [opt_single_val(j),opt_single_strategy(j,:)] = direct_chance_5time(k);
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
value = -quad_form(naive_trading,negative_P0) - dot(negative_q0,naive_trading);
plot(1.-[0:0.05:0.4,0.4:0.01:0.68],opt_single_val,'b');
yline(value,'r');
xlabel('required minimum capital ratio')
ylabel('expected revenue')
title('single asset')
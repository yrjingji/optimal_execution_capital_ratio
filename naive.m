function value_naive = naive(sigma)
    initial_price = [17,17];
    tau = 1;
    l = 1200;
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
    value_naive = -quad_form(naive_trading,negative_P0) - dot(negative_q0,naive_trading);
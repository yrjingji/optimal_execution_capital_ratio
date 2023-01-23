function [opt_val,opt_mean,opt_var,sol_first,sol_second] = general_model_portfolio_5time_meanvar_noconstraint(sigma,mu)
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
    A_obj = c1*sigma;
    B_obj = 0.5*c2*sigma;
    negative_P0 =[A_obj B_obj B_obj B_obj B_obj; B_obj A_obj B_obj B_obj B_obj;...
        B_obj B_obj A_obj B_obj B_obj;B_obj B_obj B_obj A_obj B_obj;B_obj B_obj B_obj B_obj A_obj];
    negative_q0 = -[initial_price(1);initial_price(2);initial_price(1);initial_price(2);initial_price(1);initial_price(2);...
        initial_price(1);initial_price(2);initial_price(1);initial_price(2)];
     U1 = chol(sigma);
    cvx_begin quiet
       variable s1(5)
       variable s2(5);
            %need modification
       trans_s = [s1(1),s2(1),s1(2),s2(2),s1(3),s2(3),s1(4),s2(4),s1(5),s2(5)];
            
             % for variance calculation
       q1 = sqrt(tau)*U1'*(S_0-[s1(1);s2(1)]);
       q2 = sqrt(tau)*U1'*(S_0-[s1(1);s2(1)] -[s1(2);s2(2)]);
       q3 = sqrt(tau)*U1'*([s1(4);s2(4)]+[s1(5);s2(5)]);
       q4 = sqrt(tau)*U1'*([s1(5);s2(5)]);
       var_portfolio = dot(q1,q1)+dot(q2,q2)+dot(q3,q3)+dot(q4,q4);
       minimize(quad_form(trans_s,negative_P0) + dot(negative_q0,trans_s)+mu*var_portfolio)
       subject to
            %fully invest:
           sum(s1) -S_0(1) == 0;
           sum(s2) - S_0(2) == 0;
           for i=1:4
               sum(s1(1:i)) <= X_0(1);
               sum(s2(1:i)) <= X_0(2);
            end
    cvx_end
     
    opt_val= -cvx_optval;
    opt_mean = -quad_form(trans_s,negative_P0)- dot(negative_q0,trans_s);
    opt_var = var_portfolio;
    sol_first = s1;
    sol_second = s2';
end

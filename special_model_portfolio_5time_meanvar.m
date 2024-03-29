function [opt_val,opt_mean,opt_var,sol_first,sol_second,prob_split_product] = special_model_portfolio_5time_meanvar(k,sigma,mu)
    %the decision variables need to change the order to (s1,u1,s2,u2,s3,u3,s4,u4,s5,u5)
    k1 = k(1);
    k2 = k(2);
    initial_price = [16;16];
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
   
    A = [(c1-k1*c2)*sigma(1,1), (c1-c2*(k1+k2)/2)*sigma(1,2); (c1-c2*(k1+k2)/2)*sigma(1,2), (c1-k2*c2)*sigma(2,2)];
    B = [c2*(0.5-k1)*sigma(1,1), c2*(1-k1-k2)/2*sigma(1,2);c2*(1-k1-k2)/2*sigma(1,2), c2*(0.5-k2)*sigma(2,2)];
   %first step capital ratio constraint
    negative_A1 = A;
    negative_b1 = [(k(1)-1)*initial_price(1)+k(1)*X_0(1)*gamma(1,1)+k(2)*X_0(2)*gamma(1,2),...
        (k(2)-1)*initial_price(2)+k(2)*X_0(2)*gamma(2,2) + k(1)*X_0(1)*gamma(2,1)];
    negative_c1 = l - k(1)*X_0(1)*initial_price(1) - k(2)*X_0(2)*initial_price(2);
    %first step first chance constraint's coefficients
    D1(1,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0];
    D1(2,:) = sqrt(sigma(2,2)*tau)*[0,-k(2)];
    e1 = [sqrt(sigma(1,1)*tau)*k(1)*X_0(1); sqrt(sigma(2,2)*tau)*k(2)*X_0(2)];
    %second step capital ratio constraint
    negative_A2 = [A,B;B,A];
    negative_b2 = [(k(1)-1)*initial_price(1)+k(1)*X_0(1)*gamma(1,1)+k(2)*X_0(2)*gamma(1,2),...
        (k(2)-1)*initial_price(2)+k(2)*X_0(2)*gamma(2,2) + k(1)*X_0(1)*gamma(2,1),....
        (k(1)-1)*initial_price(1)+k(1)*X_0(1)*gamma(1,1)+k(2)*X_0(2)*gamma(1,2),...
        (k(2)-1)*initial_price(2)+k(2)*X_0(2)*gamma(2,2) +  k(1)*X_0(1)*gamma(2,1)];
    negative_c2 = l - k(1)*X_0(1)*initial_price(1) - k(2)*X_0(2)*initial_price(2);
    %second step first chance constraint's coefficients
    D2(1,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0,1-k(1),0];
    D2(2,:) = sqrt(sigma(2,2)*tau)*[0,-k(2), 0,1-k(2)];
    D2(3,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0,-k(1),0];
    D2(4,:) = sqrt(sigma(2,2)*tau)*[0,-k(2), 0,-k(2)];
    e2 = [sqrt(sigma(1,1)*tau)*k(1)*X_0(1); sqrt(sigma(2,2)*tau)*k(2)*X_0(2);...
        sqrt(sigma(1,1)*tau)*k(1)*X_0(1); sqrt(sigma(2,2)*tau)*k(2)*X_0(2)];
    %third step capital ratio constraint
    negative_A3 = [A,B,B;B,A,B;B,B,A];
    negative_b3 = [(k(1)-1)*initial_price(1)+k(1)*X_0(1)*gamma(1,1)+k(2)*X_0(2)*gamma(1,2),...
        (k(2)-1)*initial_price(2)+k(2)*X_0(2)*gamma(2,2) + k(1)*X_0(1)*gamma(2,1),....
        (k(1)-1)*initial_price(1)+k(1)*X_0(1)*gamma(1,1)+k(2)*X_0(2)*gamma(1,2),...
        (k(2)-1)*initial_price(2)+k(2)*X_0(2)*gamma(2,2) +  k(1)*X_0(1)*gamma(2,1),...
        (k(1)-1)*initial_price(1)+k(1)*X_0(1)*gamma(1,1)+k(2)*X_0(2)*gamma(1,2),...
        (k(2)-1)*initial_price(2)+k(2)*X_0(2)*gamma(2,2) +  k(1)*X_0(1)*gamma(2,1)];
    negative_c3 = l - k(1)*X_0(1)*initial_price(1) - k(2)*X_0(2)*initial_price(2);
    %third step first chance constraint's coefficients
    D3(1,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0,1-k(1),0,1-k(1),0];
    D3(2,:) = sqrt(sigma(2,2)*tau)*[0,-k(2), 0,1-k(2),0,1-k(2)];
    D3(3,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0,-k(1),0,1-k(1),0];
    D3(4,:) = sqrt(sigma(2,2)*tau)*[0,-k(2), 0,-k(2),0,1-k(2)];
    D3(5,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0,-k(1),0,-k(1),0];
    D3(6,:) = sqrt(sigma(2,2)*tau)*[0,-k(2), 0,-k(2),0,-k(2)];
    e3 = [sqrt(sigma(1,1)*tau)*k(1)*X_0(1); sqrt(sigma(2,2)*tau)*k(2)*X_0(2);...
        sqrt(sigma(1,1)*tau)*k(1)*X_0(1); sqrt(sigma(2,2)*tau)*k(2)*X_0(2);...
        sqrt(sigma(1,1)*tau)*k(1)*X_0(1); sqrt(sigma(2,2)*tau)*k(2)*X_0(2)];
    
   %fourth step capital ratio constraint
    negative_A4 = [A,B,B,B;B,A,B,B;B,B,A,B;B,B,B,A];
    
    negative_b4 = [ negative_b3,[(k(1)-1)*initial_price(1)+k(1)*X_0(1)*gamma(1,1)+k(2)*X_0(2)*gamma(1,2),....
        (k(2)-1)*initial_price(2)+k(2)*X_0(2)*gamma(2,2) + k(1)*X_0(1)*gamma(2,1)]];
    negative_c4 = l - k(1)*X_0(1)*initial_price(1) - k(2)*X_0(2)*initial_price(2);
    %fourth step first chance constraint's coefficients
    D4(1,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0,1-k(1),0,1-k(1),0,1-k(1),0];
    D4(2,:) = sqrt(sigma(2,2)*tau)*[0,-k(2), 0,1-k(2),0,1-k(2),0,1-k(2)];
    D4(3,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0,-k(1),0,1-k(1),0,1-k(1),0];
    D4(4,:) = sqrt(sigma(2,2)*tau)*[0,-k(2), 0,-k(2),0,1-k(2),0,1-k(2)];
    D4(5,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0,-k(1),0,-k(1),0,1-k(1),0];
    D4(6,:) = sqrt(sigma(2,2)*tau)*[0,-k(2), 0,-k(2),0,-k(2),0,1-k(2)];
    D4(7,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0,-k(1),0,-k(1),0,-k(1),0];
    D4(8,:) = sqrt(sigma(2,2)*tau)*[0,-k(2), 0,-k(2),0,-k(2),0,-k(2)];
    
    e4 = [e3;[sqrt(sigma(1,1)*tau)*k(1)*X_0(1); sqrt(sigma(2,2)*tau)*k(2)*X_0(2)]];
    
    %fifth step capital ratio constraint
    negative_A5 = [A,B,B,B,B;B,A,B,B,B;B,B,A,B,B;B,B,B,A,B;B,B,B,B,A ];
    negative_b5 = [ negative_b4,[(k(1)-1)*initial_price(1)+k(1)*X_0(1)*gamma(1,1)+k(2)*X_0(2)*gamma(1,2),....
        (k(2)-1)*initial_price(2)+k(2)*X_0(2)*gamma(2,2) + k(1)*X_0(1)*gamma(2,1)]];
    negative_c5 = l - k(1)*X_0(1)*initial_price(1) - k(2)*X_0(2)*initial_price(2);
    %fifth step first chance constraint's coefficients
    D5(1,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0,1-k(1),0,1-k(1),0,1-k(1),0,1-k(1),0];
    D5(2,:) = sqrt(sigma(2,2)*tau)*[0,-k(2), 0,1-k(2),0,1-k(2),0,1-k(2)0,1-k(2)];
    D5(3,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0,-k(1),0,1-k(1),0,1-k(1),0,1-k(1),0];
    D5(4,:) = sqrt(sigma(2,2)*tau)*[0,-k(2), 0,-k(2),0,1-k(2),0,1-k(2),0,1-k(2)];
    D5(5,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0,-k(1),0,-k(1),0,1-k(1),0,1-k(1),0];
    D5(6,:) = sqrt(sigma(2,2)*tau)*[0,-k(2), 0,-k(2),0,-k(2),0,1-k(2),0,1-k(2)];
    D5(7,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0,-k(1),0,-k(1),0,-k(1),0,1-k(1),0];
    D5(8,:) = sqrt(sigma(2,2)*tau)*[0,-k(2), 0,-k(2),0,-k(2),0,-k(2),0,1-k(2)];
    D5(9,:) = sqrt(sigma(1,1)*tau)*[-k(1), 0,-k(1),0,-k(1),0,-k(1),0,-k(1),0];
    D5(10,:)= sqrt(sigma(2,2)*tau)*[0,-k(2), 0,-k(2),0,-k(2),0,-k(2),0,-k(2)];
    e5 = [e4;[sqrt(sigma(1,1)*tau)*k(1)*X_0(1); sqrt(sigma(2,2)*tau)*k(2)*X_0(2)]];

    %cholesky decompostion
    cov_1 = sigma;
    U1 = chol(cov_1);
    cov_2 = [sigma zeros(2,2);zeros(2,2) sigma];
    U2 = chol(cov_2);
    cov_3 = [sigma zeros(2,2) zeros(2,2);zeros(2,2) sigma  zeros(2,2);zeros(2,2) zeros(2,2) sigma];
    U3 = chol(cov_3);
    cov_4 = [sigma zeros(2,2) zeros(2,2) zeros(2,2);zeros(2,2) sigma  zeros(2,2) zeros(2,2);...
        zeros(2,2) zeros(2,2) sigma zeros(2,2); zeros(2,2) zeros(2,2) zeros(2,2) sigma];
    U4 = chol(cov_4);
    cov_5 = [sigma zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2);zeros(2,2) sigma  zeros(2,2) zeros(2,2) zeros(2,2);...
        zeros(2,2) zeros(2,2) sigma zeros(2,2) zeros(2,2); zeros(2,2) zeros(2,2) zeros(2,2) sigma zeros(2,2);
        zeros(2,2) zeros(2,2) zeros(2,2) zeros(2,2) sigma];
    U5 = chol(cov_5);
    for j=1:200
        alpha_1(j,:) = drchrnd([1,1,1],1)*0.05;
        alpha_2(j,:) = drchrnd([1,1,1],1)*0.05;
        alpha_3(j,:) = drchrnd([1,1,1],1)*0.05;
        alpha_4(j,:) = drchrnd([1,1,1],1)*0.05;
        alpha_5(j,:) = drchrnd([1,1,1],1)*0.05;
    end
    for j=1:200

        cvx_begin quiet
            variable s1(5)
            variable s2(5);
            variable b(5);
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
            %first step's first constriant
            trans_s_one = [s1(1),s2(1)];
            quad_form(trans_s_one,negative_A1)+dot(negative_b1,trans_s_one)+negative_c1+norminv(1-alpha_1(j,1))*norm(U1 * D1 *trans_s_one'+U1*e1)<=0;
            %first step second chance constraint's coefficients
            gamma(1,1)*s1(1) + gamma(2,1)*s2(1) - initial_price(1) + sqrt(sigma(1,1)*tau)*norminv(1-alpha_1(j,2)) <=0;
            gamma(1,2)*s1(1)+gamma(2,2)*s2(1) - initial_price(2) + sqrt(sigma(2,2)*tau)*norminv(1-alpha_1(j,3)) <=0;
            %second step's first constriant
            trans_s_two = [s1(1),s2(1),s1(2),s2(2)];
            quad_form(trans_s_two,negative_A2)+dot(negative_b2,trans_s_two)+negative_c2+norminv(1-alpha_2(j,1))*norm(U2 * D2 *trans_s_two'+U2*e2)<=0;
            %second step second chance constraint's coefficients
            gamma(1,1)*sum([s1(1),s1(2)])+gamma(2,1)*sum([s2(1),s2(2)]) - initial_price(1) + sqrt(sigma(1,1)*2*tau)*norminv(1-alpha_2(j,2)) <=0;
            gamma(1,2)*sum([s1(1),s1(2)])+gamma(2,2)*sum([s2(1),s2(2)]) - initial_price(2) + sqrt(sigma(2,2)*2*tau)*norminv(1-alpha_2(j,3)) <=0;
             %third step's first constriant
            trans_s_third = [s1(1),s2(1),s1(2),s2(2),s1(3),s2(3)];
            quad_form(trans_s_third,negative_A3)+dot(negative_b3,trans_s_third)+negative_c3+norminv(1-alpha_3(j,1))*norm(U3 * D3 *trans_s_third'+U3*e3)<=0;
            %third step second chance constraint's coefficients
            gamma(1,1)*sum([s1(1),s1(2),s1(3)])+gamma(2,1)*sum([s2(1),s2(2),s2(3)]) - initial_price(1) + sqrt(sigma(1,1)*3*tau)*norminv(1-alpha_3(j,2)) <=0;
            gamma(1,2)*sum([s1(1),s1(2),s1(3)])+gamma(2,2)*sum([s2(1),s2(2),s2(3)]) - initial_price(2) + sqrt(sigma(2,2)*3*tau)*norminv(1-alpha_3(j,3)) <=0;
             %fourth step's first constriant
            trans_s_fourth = [s1(1),s2(1),s1(2),s2(2),s1(3),s2(3),s1(4),s2(4)];
            quad_form(trans_s_fourth,negative_A4)+dot(negative_b4,trans_s_fourth)+negative_c4+norminv(1-alpha_4(j,1))*norm(U4 * D4 *trans_s_fourth'+U4*e4)<=0;
            %fourth step second chance constraint's coefficients
            gamma(1,1)*sum([s1(1),s1(2),s1(3),s1(4)])+gamma(2,1)*sum([s2(1),s2(2),s2(3),s2(4)]) - initial_price(1) + sqrt(sigma(1,1)*4*tau)*norminv(1-alpha_4(j,2)) <=0;
            gamma(1,2)*sum([s1(1),s1(2),s1(3),s1(4)])+gamma(2,2)*sum([s2(1),s2(2),s2(3),s2(4)]) - initial_price(2) + sqrt(sigma(2,2)*4*tau)*norminv(1-alpha_4(j,3)) <=0;
            %fifth step's first constraint
            trans_s_fifth= [s1(1),s2(1),s1(2),s2(2),s1(3),s2(3),s1(4),s2(4),s1(5),s2(5)];
            quad_form(trans_s_fifth,negative_A5)+dot(negative_b5,trans_s_fifth)+negative_c5+norminv(1-alpha_5(j,1))*norm(U5 * D5 *trans_s_fifth'+U5*e5)<=0;

            %fifth step's second chance constraint
            gamma(1,1)*sum(s1)+gamma(2,1)*sum(s2) - initial_price(1) + sqrt(sigma(1,1)*5*tau)*norminv(1-alpha_5(j,2)) <=0;
            gamma(1,2)*sum(s1)+gamma(2,2)*sum(s2) - initial_price(2) + sqrt(sigma(2,2)*5*tau)*norminv(1-alpha_5(j,3)) <=0;
           %fully invest:
            sum(s1) -S_0(1) == 0;
            sum(s2) - S_0(2) == 0;
            for i=1:4
                sum(s1(1:i)) <= X_0(1);
                sum(s2(1:i)) <= X_0(2);
            end
            sum(b) == 1;
            s1 == S_0(1)*b;
            s2 == S_0(2)*b;
        cvx_end
        val(j) = -cvx_optval;
        mean_port(j) = -quad_form(trans_s,negative_P0)- dot(negative_q0,trans_s);
        var_port(j) = var_portfolio;
        sol1(j,:)= s1';
        sol2(j,:) = s2';
        end
    [opt_val,index] = max(val);
    opt_mean = mean_port(index);
    opt_var = var_port(index);
    sol_first = sol1(index,:);
    sol_second = sol2(index,:);
   prob_split_product = [alpha_1(index,1),alpha_2(index,1),alpha_3(index,1),alpha_4(index,1)];
end

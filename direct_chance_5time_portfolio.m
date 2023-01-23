function [opt_val,mean,variance,sol_first,sol_second,sol_third]=direct_chance_5time_portfolio(k,mu)
    %no cross asset impact, only diagonal matrix in price impact matrix.
    initial_price = [17,17,17];
    sigma = [0.7 0.1 -0.1; 0.1 0.7 0;-0.1 0 0.7];
    
    tau = 1;
    l = [1200,1200,1200];
    X_0 = 305*ones(1,3);
    S_0 = 250*ones(1,3);
    beta = diag(.09*ones(1,3));
    gamma = diag(.08*ones(1,3));    
    for i=1:3
        %objective function
        negative_obj_P{i} = 0.5*gamma(i,i)*ones(5,5) + diag(ones(1,5)*(beta(i,i) - 0.5*gamma(i,i)));
        negative_obj_q{i} = -initial_price(i)*ones(1,5);
        % first constraint
        negative_A1{i} = beta(i,i)-k(i)*gamma(i,i);
        negative_b1{i} = (k(i)-1)*initial_price(i) + k(i)*X_0(i)*gamma(i,i);
        negative_c1{i} = l(i)- k(i)*X_0(i)*initial_price(i);

        % second constraint:
        negative_A2{i} = ones(2,2)*(0.5-k(i))*gamma(i,i) + diag([beta(i,i)- k(i)*gamma(i,i) - (0.5-k(i))*gamma(i,i),beta(i,i)- k(i)*gamma(i,i) - (0.5-k(i))*gamma(i,i)]);
        negative_b2{i} = [(k(i)-1)*initial_price(i) + k(i)*X_0(i)*gamma(i,i),(k(i)-1)*initial_price(i) + k(i)*X_0(i)*gamma(i,i)];
        negative_c2{i} = l(i)- k(i)*X_0(i)*initial_price(i);
       
        % third constraint:
        negative_A3{i} = ones(3,3)*(0.5-k(i))*gamma(i,i) + diag((beta(i,i)- k(i)*gamma(i,i) - (0.5-k(i))*gamma(i,i))*ones(1,3));
        negative_b3{i} = ((k(i)-1)*initial_price(i) + k(i)*X_0(i)*gamma(i,i))*ones(1,3);
        negative_c3{i} = l(i)- k(i)*X_0(i)*initial_price(i);
       
        % fourth consitraint:
        negative_A4{i}= ones(4,4)*(0.5-k(i))*gamma(i,i) + diag((beta(i,i)- k(i)*gamma(i,i) - (0.5-k(i))*gamma(i,i))*ones(1,4));
        negative_b4{i} = ((k(i)-1)*initial_price(i) + k(i)*X_0(i)*gamma(i,i))*ones(1,4);
        negative_c4{i} = l(i)- k(i)*X_0(i)*initial_price(i);
       
    end
    %coefficients of objective function
    O1(1,1:5) = [0 1 1 1 1];
    O1(2,6:10) = [0 1 1 1 1];
    O1(3,11:15) = [0 1 1 1 1];
     O2(1,1:5) = [0 0 1 1 1];
    O2(2,6:10) = [0 0 1 1 1];
    O2(3,11:15) = [0 0 1 1 1];
     O3(1,1:5) = [0 0 0 1 1];
    O3(2,6:10) = [0 0 0 1 1];
    O3(3,11:15) = [0 0 0 1 1];
     O4(1,1:5) = [0 0 0 0 1];
    O4(2,6:10) = [0 0 0 0 1];
    O4(3,11:15) = [0 0 0 0 1];
    %coefficients:
    %step 1
    D1 = diag([-sigma(1,1)*sqrt(tau)*k(1),-sigma(2,2)*sqrt(tau)*k(2),-sigma(3,3)*sqrt(tau)*k(3)]);
    e1 = [sigma(1,1)*sqrt(tau)*k(1)*X_0(1);sigma(2,2)*sqrt(tau)*k(2)*X_0(2);sigma(3,3)*sqrt(tau)*k(3)*X_0(3)];
    %step 2
    i = 1;
    for j=1:2:5
        D2(j:j+1,j:j+1) = sigma(i,i)*sqrt(tau)*[-k(i), 1-k(i),;-k(i), -k(i)];
        i= i+1;
    end
    e2 = [sigma(1,1)*sqrt(tau)*k(1)*X_0(1)*ones(2,1);sigma(2,2)*sqrt(tau)*k(2)*X_0(2)*ones(2,1);sigma(3,3)*sqrt(tau)*k(3)*X_0(3)*ones(2,1)];
    % step 3
    i=1;
    for j=1:3:7
        D3(j:j+2,j:j+2) = sigma(i,i)*sqrt(tau)*[-k(i), 1-k(i),1-k(i);-k(i), -k(i),1-k(i);-k(i),-k(i),-k(i)];
        i = i+1;
    end
    e3 = [sigma(1,1)*sqrt(tau)*k(1)*X_0(1)*ones(3,1);sigma(2,2)*sqrt(tau)*k(2)*X_0(2)*ones(3,1);sigma(3,3)*sqrt(tau)*k(3)*X_0(3)*ones(3,1)];
    % step 4
    i=1;
    for j=1:4:9
        D4(j:j+3,j:j+3) = sigma(i,i)*sqrt(tau)*[-k(i), 1-k(i),1-k(i),1-k(i);-k(i), -k(i),1-k(i),1-k(i);-k(i),-k(i),-k(i),1-k(i);-k(i),-k(i),-k(i),-k(i)];
        i = i+1;
    end
    e4 = [sigma(1,1)*sqrt(tau)*k(1)*X_0(1)*ones(4,1);sigma(2,2)*sqrt(tau)*k(2)*X_0(2)*ones(4,1);sigma(3,3)*sqrt(tau)*k(3)*X_0(3)*ones(4,1)];
    %covariance matrix needs updated, as well as the Cholesky
    %decomposition.
    cov_1 = sigma;
    upper_2 = diag([sigma(1,1)*ones(1,2), sigma(2,2)*ones(1,2),sigma(3,3)*ones(1,2)])+...
        diag([sigma(1,2)*ones(1,2),sigma(2,3)*ones(1,2)],2)+ diag([sigma(1,3)*ones(1,2)],4);
    cov_2 = triu(upper_2,1)'+upper_2;
    upper_3 = diag([sigma(1,1)*ones(1,3), sigma(2,2)*ones(1,3),sigma(3,3)*ones(1,3)])+...
        diag([sigma(1,2)*ones(1,3),sigma(2,3)*ones(1,3)],3)+ diag([sigma(1,3)*ones(1,3)],6);
    cov_3 = triu(upper_3,1)'+upper_3;
    upper_4 = diag([sigma(1,1)*ones(1,4), sigma(2,2)*ones(1,4),sigma(3,3)*ones(1,4)])+...
        diag([sigma(1,2)*ones(1,4),sigma(2,3)*ones(1,4)],4)+ diag([sigma(1,3)*ones(1,4)],8);
    
    cov_4 = triu(upper_4,1)'+upper_4;
    %cholesky decompostion
    U_1 = chol(cov_1);
    U_2 = chol(cov_2);
    U_3 = chol(cov_3);
    U_4 = chol(cov_4);
    rng('default');
    rng(1);
    for j=1:200
        
        alpha_1 = drchrnd([1,1,1,1],1)*0.05;
        alpha_2 = drchrnd([1,1,1,1],1)*0.05;
        alpha_3 = drchrnd([1,1,1,1],1)*0.05;
        alpha_4 = drchrnd([1,1,1,1],1)*0.05;
        cvx_begin quiet
            variable s1(5)
            variable s2(5)
            variable s3(5)
            d1 = O1*[s1;s2;s3];
            d2 =O2*[s1;s2;s3];
            d3 = O3*[s1;s2;s3];
            d4 = O4*[s1;s2;s3];
            minimize(quad_form(s1,negative_obj_P{1}) + quad_form(s2,negative_obj_P{2})+quad_form(s3,negative_obj_P{3}) + dot(negative_obj_q{1},s1)+dot(negative_obj_q{2},s2)+dot(negative_obj_q{3},s3)...
                +mu*(d1'*sigma*d1 + d2'*sigma*d2 + d3'*sigma*d3 + d4'*sigma*d4))
            subject to
            %first step
            quad_form(s1(1),negative_A1{1})+quad_form(s2(1),negative_A1{2})+ quad_form(s3(1),negative_A1{3})+ dot(negative_b1{1},s1(1))+ ...
                dot(negative_b1{2},s2(1)) + dot(negative_b1{3},s3(1))+negative_c1{1}+...
                negative_c1{2}+negative_c1{3}+ norminv(1-alpha_1(1))*norm(U_1 * D1 *([s1(1);s2(1);s3(1)])+U_1*e1)<=0;
            gamma(1,1)*s1(1) - initial_price(1) + sigma(1,1)*sqrt(tau)*norminv(1-alpha_1(2)) <=0;
            gamma(2,2)*s2(1) - initial_price(2) + sigma(2,2)*sqrt(tau)*norminv(1-alpha_1(3)) <=0;
            gamma(3,3)*s3(1) - initial_price(3) + sigma(3,3)*sqrt(tau)*norminv(1-alpha_1(4)) <=0;
            %second step
            quad_form([s1(1),s1(2)],negative_A2{1})+quad_form([s2(1),s2(2)],negative_A2{2})+ quad_form([s3(1),s3(2)],negative_A2{3})+ dot(negative_b2{1},[s1(1),s1(2)])+ ...
                dot(negative_b2{2},[s2(1),s2(2)]) + dot(negative_b2{3},[s3(1),s3(2)])+negative_c2{1}+...
                negative_c2{2}+negative_c2{3}+ norminv(1-alpha_2(1))*norm(U_2 * D2 *([s1(1);s1(2);s2(1);s2(2);s3(1);s3(2)])+U_2*e2)<=0;
            gamma(1,1)*sum([s1(1),s1(2)]) - initial_price(1) + sigma(1,1)*sqrt(2*tau)*norminv(1-alpha_2(2)) <=0;
            gamma(2,2)*sum([s2(1),s2(2)]) - initial_price(2) + sigma(2,2)*sqrt(2*tau)*norminv(1-alpha_2(3)) <=0;
            gamma(3,3)*sum([s3(1),s3(2)]) - initial_price(3) + sigma(3,3)*sqrt(2*tau)*norminv(1-alpha_2(4)) <=0;
            %third step
            
            quad_form([s1(1),s1(2),s1(3)],negative_A3{1})+quad_form([s2(1),s2(2),s2(3)],negative_A3{2})+ quad_form([s3(1),s3(2),s3(3)],negative_A3{3})+ dot(negative_b3{1},[s1(1),s1(2),s1(3)])+ ...
                dot(negative_b3{2},[s2(1),s2(2),s2(3)]) + dot(negative_b3{3},[s3(1),s3(2),s3(3)])+negative_c3{1}+...
                negative_c3{2}+negative_c3{3}+ norminv(1-alpha_3(1))*norm(U_3 * D3 *([s1(1);s1(2);s1(3);s2(1);s2(2);s2(3);s3(1);s3(2);s3(3)])+U_3*e3)<=0;
            gamma(1,1)*sum([s1(1),s1(2),s1(3)]) - initial_price(1) + sigma(1,1)*sqrt(3*tau)*norminv(1-alpha_3(2)) <=0;
            gamma(2,2)*sum([s2(1),s2(2),s2(3)]) - initial_price(2) + sigma(2,2)*sqrt(3*tau)*norminv(1-alpha_3(3)) <=0;
            gamma(3,3)*sum([s3(1),s3(2),s3(3)]) - initial_price(3) + sigma(3,3)*sqrt(3*tau)*norminv(1-alpha_3(4)) <=0;
            %fourth step
            quad_form([s1(1),s1(2),s1(3),s1(4)],negative_A4{1})+quad_form([s2(1),s2(2),s2(3),s2(4)],negative_A4{2})+ quad_form([s3(1),s3(2),s3(3),s3(4)],negative_A4{3})+ dot(negative_b4{1},[s1(1),s1(2),s1(3),s1(4)])+ ...
                dot(negative_b4{2},[s2(1),s2(2),s2(3),s2(4)]) + dot(negative_b4{3},[s3(1),s3(2),s3(3),s3(4)])+negative_c4{1}+...
                negative_c4{2}+negative_c4{3}+ norminv(1-alpha_4(1))*norm(U_4 * D4 *([s1(1);s1(2);s1(3);s1(4);s2(1);s2(2);s2(3);s2(4);s3(1);s3(2);s3(3);s3(4)])+U_4*e4)<=0;
            gamma(1,1)*sum([s1(1),s1(2),s1(3),s1(4)]) - initial_price(1) + sigma(1,1)*sqrt(4*tau)*norminv(1-alpha_4(2)) <=0;
            gamma(2,2)*sum([s2(1),s2(2),s2(3),s2(4)]) - initial_price(2) + sigma(2,2)*sqrt(4*tau)*norminv(1-alpha_4(3)) <=0;
            gamma(3,3)*sum([s3(1),s3(2),s3(3),s3(4)]) - initial_price(3) + sigma(3,3)*sqrt(4*tau)*norminv(1-alpha_4(4)) <=0;
            
            %fully invest:
            sum(s1) -S_0(1) == 0;
            sum(s2) - S_0(2) == 0;
            sum(s3) - S_0(3) ==0;
            for i=1:5
                sum(s1(1:i)) <= X_0(1);
                sum(s2(1:i)) <= X_0(2);
                sum(s3(1:i)) <= X_0(3);
            end
        cvx_end
        val(j) = -cvx_optval;
        mean_val(j) = quad_form(s1,negative_obj_P{1}) + quad_form(s2,negative_obj_P{2})+quad_form(s3,negative_obj_P{3}) + dot(negative_obj_q{1},s1)+dot(negative_obj_q{2},s2)+dot(negative_obj_q{3},s3);
        variance_val(j) = d1'*sigma*d1 + d2'*sigma*d2 + d3'*sigma*d3 + d4'*sigma*d4;
        sol_1(j,:)= s1';
        sol_2(j,:) = s2';
        sol_3(j,:) = s3';
        end
      [opt_val,index] = max(val);
      mean = mean_val(index);
      variance = variance_val(index);
      sol_first = sol_1(index,:);
      sol_second = sol_2(index,:);
      sol_third = sol_3(index,:);
end



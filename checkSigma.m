
parfor j=1:length(all_combinations)
    if all_combinations(j,7) > all_combinations(j,8) 
            [~,~,~,optSol(j,:),~] = direct_chance_5time_mean_var(all_combinations(j,:),0.5,0);

    end
end

%record the probability split to the product one 
function [opt_val,mean_val,variance_val,opt_sol,prob_split_product]=direct_chance_5time_mean_var(params,k,mu)
    sigma = params[1];
    initial_price = params[2];
    %standard deviation
    tau = params[3];
    l = params[4];
    X_0 = params[5];
    S_0 = params[6];
    beta = params[7];
    gamma = params[8];
    negative_P0 = 0.5*gamma*ones(5,5) + diag(ones(1,5)*(beta - 0.5*gamma));
    negative_q0 = -initial_price*ones(1,5);
    %first step:
    negative_A1 = beta-k*gamma;
    negative_b1 = (k-1)*initial_price + k*X_0*gamma;
    negative_c1 = l- k*X_0*initial_price;
    %first step coefficients:
    D1 = -k;
    e1= k*X_0;
    %second step:
    negative_A2 = ones(2,2)*(0.5-k)*gamma + diag([beta- k*gamma - (0.5-k)*gamma,beta- k*gamma - (0.5-k)*gamma]);
    negative_b2 = [(k-1)*initial_price + k*X_0*gamma,(k-1)*initial_price + k*X_0*gamma];
    negative_c2 = l - k*X_0*initial_price;
    %second step:
    D2 = [-k, 1-k;-k, -k];
    e2 = [k*X_0;k*X_0];
     %third step:
    negative_A3 = ones(3,3)*(0.5-k)*gamma + diag((beta- k*gamma - (0.5-k)*gamma)*ones(1,3));
    negative_b3 = ((k-1)*initial_price + k*X_0*gamma)*ones(1,3);
    negative_c3 = l - k*X_0*initial_price;
    %third step:
    D3 = [-k, 1-k,1-k;-k, -k,1-k;-k,-k,-k];
    e3 = [k*X_0;k*X_0;k*X_0];
    
     %fourth step:
    negative_A4 = ones(4,4)*(0.5-k)*gamma + diag((beta- k*gamma - (0.5-k)*gamma)*ones(1,4));
    negative_b4 = ((k-1)*initial_price + k*X_0*gamma)*ones(1,4);
    negative_c4 = l - k*X_0*initial_price;
    %fourth step:
    D4 = [-k, 1-k,1-k,1-k;-k, -k,1-k,1-k;-k,-k,-k,1-k;-k,-k,-k,-k];
    e4 = k*X_0*ones(4,1);
    %fifth step:
    negative_A5 = ones(5,5)*(0.5-k)*gamma + diag((beta- k*gamma - (0.5-k)*gamma)*ones(1,5));
    negative_b5 = ((k-1)*initial_price + k*X_0*gamma)*ones(1,5);
    negative_c5 = l - k*X_0*initial_price;
    %fifth step:
    D5 = [-k, 1-k,1-k,1-k, 1-k;-k, -k,1-k,1-k, 1-k;-k,-k,-k,1-k, 1-k;-k,-k,-k,-k, 1-k; -k*ones(1,5)];
    e5 = k*X_0*ones(5,1);

    %coefficients for variance in objective.
    Q(1,:) = sqrt(tau)*sigma*[0 1 1 1 1];
    Q(2,:) = sqrt(tau)*sigma*[0 0 1 1 1];
    Q(3,:) = sqrt(tau)*sigma*[0 0 0 1 1];
    Q(4,:) = sqrt(tau)*sigma*[0 0 0 0 1];
    for j=1:500
        alpha_1(j,:) = drchrnd([1,1],1)*0.05;
        alpha_2(j,:) = drchrnd([1,1],1)*0.05;
        alpha_3(j,:) = drchrnd([1,1],1)*0.05;
        alpha_4(j,:) = drchrnd([1,1],1)*0.05;
        alpha_5(j,:) = drchrnd([1,1],1)*0.05;
    end
    for j=1:500
        cvx_begin quiet
            variable s(5)
            var_coeff = Q*s;
            minimize(quad_form(s,negative_P0) + dot(negative_q0,s)+mu*(var_coeff'*var_coeff))
            subject to
        %first step
                quad_form(s(1),negative_A1) + dot(negative_b1,s(1)) + negative_c1 + sigma*sqrt(tau)*norminv(1-alpha_1(j,1))*norm(D1*s(1) + e1)<=0;
                gamma*s(1) - initial_price + sigma*sqrt(tau)*norminv(1-alpha_1(j,2)) <=0;
        %second step
                quad_form([s(1);s(2)],negative_A2) + dot(negative_b2,[s(1);s(2)]) + negative_c2 + sigma*sqrt(tau)*norminv(1-alpha_2(j,1))*norm(D2*[s(1);s(2)] + e2)<=0;
                gamma*sum([s(1);s(2)]) - initial_price + sigma*sqrt(2*tau)*norminv(1-alpha_2(j,2)) <=0;
        %third step
                quad_form([s(1);s(2);s(3)],negative_A3) + dot(negative_b3,[s(1);s(2);s(3)]) + negative_c3 + sigma*sqrt(tau)*norminv(1-alpha_3(j,1))*norm(D3*[s(1);s(2);s(3)] + e3)<=0;
                gamma*sum([s(1);s(2);s(3)]) - initial_price + sigma*sqrt(3*tau)*norminv(1-alpha_3(j,2)) <=0;
        %fourth step
                quad_form([s(1);s(2);s(3);s(4)],negative_A4) + dot(negative_b4,[s(1);s(2);s(3);s(4)]) + negative_c4 + sigma*sqrt(tau)*norminv(1-alpha_4(j,1))*norm(D4*[s(1);s(2);s(3);s(4)] + e4)<=0;
                gamma*sum([s(1);s(2);s(3);s(4)]) - initial_price + sigma*sqrt(4*tau)*norminv(1-alpha_4(j,2)) <=0;
        %fifth step
                quad_form(s,negative_A5) + dot(negative_b5,s) + negative_c5 + sigma*sqrt(tau)*norminv(1-alpha_5(j,1))*norm(D5*s + e5)<=0;
                gamma*sum(s) - initial_price + sigma*sqrt(5*tau)*norminv(1-alpha_5(j,2)) <=0;
                
         %fully invest:
                sum(s) - S_0 == 0;
                for i=1:5
                    sum(s(1:i)) <= X_0;
                end
        cvx_end
        val(j) = -cvx_optval;
        mean(j) = -quad_form(s,negative_P0) - dot(negative_q0,s);
        variance(j) = var_coeff'*var_coeff;
        sol(j,1:5)= s';
    end
    [opt_val,index] = max(val);
    mean_val = mean(index);
    variance_val = variance(index);
    opt_sol = sol(index,:);
    prob_split_product = [alpha_1(index,1),alpha_2(index,1),alpha_3(index,1),alpha_4(index,1)];
end


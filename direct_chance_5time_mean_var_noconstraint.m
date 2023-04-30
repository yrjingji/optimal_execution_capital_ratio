function [opt_val,mean_val,variance_val,opt_sol]= direct_chance_5time_mean_var_noconstraint(mu)
    initial_price = 16;
    %standard deviation
    sigma = 0.77;
    tau = 1;
    l = 1200;
    X_0 = 305;
    S_0 = 250;
    beta = .095;
    gamma = .05;
    negative_P0 = 0.5*gamma*ones(5,5) + diag(ones(1,5)*(beta - 0.5*gamma));
    negative_q0 = -initial_price*ones(1,5);
    %coefficients for variance in objective.
    Q(1,:) = sqrt(tau)*sigma*[0 1 1 1 1];
    Q(2,:) = sqrt(tau)*sigma*[0 0 1 1 1];
    Q(3,:) = sqrt(tau)*sigma*[0 0 0 1 1];
    Q(4,:) = sqrt(tau)*sigma*[0 0 0 0 1];

    cvx_begin quiet
         variable s(5)
         var_coeff = Q*s;
         minimize(quad_form(s,negative_P0) + dot(negative_q0,s)+mu*(var_coeff'*var_coeff))
         subject to
        %fully invest:
            sum(s) - S_0 == 0;
            for i=1:5
                sum(s(1:i)) <= X_0;
            end
    cvx_end
    opt_val = -cvx_optval;
    mean_val = -quad_form(s,negative_P0) - dot(negative_q0,s);
    variance_val = var_coeff'*var_coeff;
    opt_sol= s';
end


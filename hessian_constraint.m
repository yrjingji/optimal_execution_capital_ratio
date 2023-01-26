function [first,second,H] = hessian_constraint(q)
    initial_price = [17,17];
    sigma = [.007, 0.00006;0.00006,.00007];
    U_2 = chol(sigma);
    tau = 1;
    l = 1200;
    X_0 = [305,305];
    S_0 = [300,300];
    beta = 0.01*sigma;
    gamma = 0.009*sigma;
    s = (S_0/2)';
    risk_weights = [8,8];
    B = diag([risk_weights(1)*gamma(1,1),risk_weights(2)*gamma(2,2)]);
    %cross impact
    B_cross_upper = diag([risk_weights(1)*gamma(1,1),risk_weights(2)*gamma(2,2)])+diag(0.5*[risk_weights(1)*gamma(2,1)+risk_weights(2)*gamma(1,2)],1);
    B_cross= triu(B_cross_upper,1)'+B_cross_upper;
    
    e_2 = [risk_weights(1)*X_0(1)*sigma(1,1);risk_weights(2)*X_0(2)*sigma(2,2)];
    P = diag([-risk_weights(1), -risk_weights(2)]);
    Q = sqrt(tau)*U_2*P;
    e = U_2*e_2;
    first = -2*B;
    norm_hessian = Q'*Q/norm(Q*s +e) - (Q'*Q*s + Q'*e)*(s'*Q'*Q + e'*Q)/(norm(Q*s+e)^3);
    second = norminv(q)*norm_hessian;
    Hessian = first+second;
    H = eig(Hessian);
end
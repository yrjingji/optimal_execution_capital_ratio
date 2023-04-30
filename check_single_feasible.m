function real_q = check_single_feasible(s,k)
    initial_price = 16;
    %standard deviation
    sigma = 0.77;
    tau = 1;
    l = 1200;
    X_0 = 305;
    S_0 = 250;
    beta = .095;
    gamma = .05;
    %check the violation of naive trading strategy for each time step
    %1,2,3,4
    s1 =s(1);
    s2 = s(1:2);
    s3 = s(1:3);
    s4 = s(1:4);
    % at the end of first period
    negative_A1 = beta-k*gamma;
    negative_b1 = (k-1)*initial_price + k*X_0*gamma;
    negative_c1 = l- k*X_0*initial_price;
    %first step coefficients:
    D1 = -k;
    e1= k*X_0;
    %first part's vector
    t1 = sigma*sqrt(tau)*(e1+ D1*s1);
    q1 = sigma*sqrt(tau)*1;
    %coordinates of the point we want to evaluate
    p_1 = quad_form(s1,negative_A1) + dot(negative_b1,s1) + negative_c1;
    p_2 = gamma*sum(s1) - initial_price;
    %calculate the probability
    mu_1 = 0;
    sigma_1 = sqrt(dot(t1,t1));
    mu_2 = 0;
    sigma_2 = sqrt(dot(q1,q1));   
    cov1 = dot(t1,q1);
    real_q(1) = 1 - normcdf(p_1,mu_1,sigma_1)- normcdf(p_2,mu_2,sigma_2)+mvncdf([p_1 p_2],[mu_1,mu_2],[sigma_1^2 cov1;cov1 sigma_2^2]+1e-6*eye(2));
    
    
    %second step:
    negative_A2 = ones(2,2)*(0.5-k)*gamma + diag([beta- k*gamma - (0.5-k)*gamma,beta- k*gamma - (0.5-k)*gamma]);
    negative_b2 = [(k-1)*initial_price + k*X_0*gamma,(k-1)*initial_price + k*X_0*gamma];
    negative_c2 = l - k*X_0*initial_price;
    %second step:
    D2 = [-k, 1-k;-k, -k];
    e2 = [k*X_0;k*X_0];
    
    %first part's vector
    t2 = sigma*sqrt(tau)*(e2+ D2*s2);
    q2 = sigma*sqrt(tau)*ones(1,2);
    %coordinates of the point we want to evaluate
    p_1 = quad_form(s2,negative_A2) + dot(negative_b2,s2) + negative_c2;
    p_2 = gamma*sum(s2) - initial_price;
    %calculate the probability
    mu_1 = 0;
    sigma_1 = sqrt(dot(t2,t2));
    mu_2 = 0;
    sigma_2 = sqrt(dot(q2,q2));   
    cov1 = dot(t2,q2);
    real_q(2) = 1 - normcdf(p_1,mu_1,sigma_1)- normcdf(p_2,mu_2,sigma_2)+mvncdf([p_1 p_2],[mu_1,mu_2],[sigma_1^2 cov1;cov1 sigma_2^2]+1e-6*eye(2));
    
    %third step:
    negative_A3 = ones(3,3)*(0.5-k)*gamma + diag((beta- k*gamma - (0.5-k)*gamma)*ones(1,3));
    negative_b3 = ((k-1)*initial_price + k*X_0*gamma)*ones(1,3);
    negative_c3 = l - k*X_0*initial_price;
    %third step:
    D3 = [-k, 1-k,1-k;-k, -k,1-k;-k,-k,-k];
    e3 = [k*X_0;k*X_0;k*X_0];
    %first part's vector
    t3 = sigma*sqrt(tau)*(e3+ D3*s3);
    q3 = sigma*sqrt(tau)*ones(1,3);
    %coordinates of the point we want to evaluate
    p_1 = quad_form(s3,negative_A3) + dot(negative_b3,s3) + negative_c3;
    p_2 = gamma*sum(s3) - initial_price;
    %calculate the probability
    mu_1 = 0;
    sigma_1 = sqrt(dot(t3,t3));
    mu_2 = 0;
    sigma_2 = sqrt(dot(q3,q3));   
    cov1 = dot(t3,q3);
    real_q(3) = 1 - normcdf(p_1,mu_1,sigma_1)- normcdf(p_2,mu_2,sigma_2)+mvncdf([p_1 p_2],[mu_1,mu_2],[sigma_1^2 cov1;cov1 sigma_2^2]+1e-6*eye(2));
     %fourth step:
    negative_A4 = ones(4,4)*(0.5-k)*gamma + diag((beta- k*gamma - (0.5-k)*gamma)*ones(1,4));
    negative_b4 = ((k-1)*initial_price + k*X_0*gamma)*ones(1,4);
    negative_c4 = l - k*X_0*initial_price;
    %fourth step:
    D4 = [-k, 1-k,1-k,1-k;-k, -k,1-k,1-k;-k,-k,-k,1-k;-k,-k,-k,-k];
    e4 = k*X_0*ones(4,1);
    %first part's vector
    t4 = sigma*sqrt(tau)*(e4+ D4*s4);
    q4 = sigma*sqrt(tau)*ones(1,4);
    %coordinates of the point we want to evaluate
    p_1 = quad_form(s4,negative_A4) + dot(negative_b4,s4) + negative_c4;
    p_2 = gamma*sum(s4) - initial_price;
    %calculate the probability
    mu_1 = 0;
    sigma_1 = sqrt(dot(t4,t4));
    mu_2 = 0;
    sigma_2 = sqrt(dot(q4,q4));   
    cov1 = dot(t4,q4);
    real_q(4) = 1 - normcdf(p_1,mu_1,sigma_1)- normcdf(p_2,mu_2,sigma_2)+mvncdf([p_1 p_2],[mu_1,mu_2],[sigma_1^2 cov1;cov1 sigma_2^2]+1e-6*eye(2));
    
    

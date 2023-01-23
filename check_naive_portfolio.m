function real_q = check_naive_portfolio(k,sigma)
%change the covariance matrix structure will make the results different a
%lot
    vol_matrix = chol(sigma);
    
    initial_price = [17;17];
    tau = 1;
    l = 1200;
    X_0 = [305;305];
    S_0 = [250;200];
    c1 = .07;
    c2 = 0.04;
    beta = c1*sigma;
    gamma = c2*sigma;
    %create holding vectors over time.

    S(1,:) = S_0(1)/5*ones(1,4);
    S(2,:) = S_0(2)/5*ones(1,4);
    X(:,1) = X_0;
    for i=2:5
        X(:,i) = X(:,i-1) - S(:,i-1);
    end
    %delete the first column
    X(:,1) = []; 
    % risk weights
    risk_weights = 8*ones(1,2);
    %at the end of first period
    t1 = sqrt(tau)*vol_matrix*diag(k)*X(:,1);
    var_1 = dot(t1,t1);
    q1 = sqrt(tau)*vol_matrix*diag(risk_weights)*X(:,1);
    var_2 = dot(q1,q1);
    cov_1 = dot(t1,q1);
    %coordinates of the point we want to evaluate
    p_1 = l-(initial_price-beta*S(:,1))'*S(:,1) - (initial_price - gamma*S(:,1))'*diag(k)*X(:,1);
    p_2 = -(initial_price - gamma*S(:,1))'*diag(risk_weights)*X(:,1);
    
    real_q(1) = 1 - normcdf(p_1,0,sqrt(var_1))- normcdf(p_2,0,sqrt(var_2))+mvncdf([p_1 p_2],[0,0],[var_1 cov_1;cov_1 var_2]+1e-6*eye(2));
    
    %at the end of second period
    t11 = sqrt(tau)*vol_matrix*diag(k)*X(:,2) + sqrt(tau)*vol_matrix*S(:,2);
    t12 = sqrt(tau)*vol_matrix*diag(k)*X(:,2);
    t1 = [t11;t12];
    var_1 = dot(t1,t1);
    q11 = sqrt(tau)*vol_matrix*diag(risk_weights)*X(:,2);
    q1 = [q11;q11];
    var_2 = dot(q1,q1);
    cov_1 = dot(t1,q1);
    %coordinates of the point we want to evaluate
    p_1 = l-(initial_price-beta*S(:,1))'*S(:,1) - (initial_price - gamma*S(:,1)-beta*S(:,2))'*S(:,2)-...
        (initial_price-gamma*(S(:,1) + S(:,2)))'*diag(k)*X(:,2);
    p_2 = -(initial_price - gamma*(S(:,1)+S(:,2)))'*diag(risk_weights)*X(:,2);
    
    real_q(2) = 1 - normcdf(p_1,0,sqrt(var_1))- normcdf(p_2,0,sqrt(var_2))+mvncdf([p_1 p_2],[0,0],[var_1 cov_1;cov_1 var_2]+1e-6*eye(2));
    
    % at the end of third step
    t11 = sqrt(tau)*vol_matrix*diag(k)*X(:,3) + sqrt(tau)*vol_matrix*(S(:,2)+S(:,3));
    t12 = sqrt(tau)*vol_matrix*diag(k)*X(:,3) + sqrt(tau)*vol_matrix*S(:,3);
    t13 = sqrt(tau)*vol_matrix*diag(k)*X(:,3);
    t1 = [t11;t12;t13];
    var_1 = dot(t1,t1);
    q11 = sqrt(tau)*vol_matrix*diag(risk_weights)*X(:,3);
    q1 = [q11;q11;q11];
    var_2 = dot(q1,q1);
    cov_1 = dot(t1,q1);
    %coordinates of the point we want to evaluate
    p_1 = l-(initial_price-beta*S(:,1))'*S(:,1) - (initial_price - gamma*S(:,1)-beta*S(:,2))'*S(:,2)-...
        (initial_price - gamma*(S(:,1)+S(:,2))-beta*S(:,3))'*S(:,3)-(initial_price-gamma*(S(:,1) + S(:,2)+S(:,3)))'*diag(k)*X(:,3);
    p_2 = -(initial_price - gamma*(S(:,1)+S(:,2)+S(:,3)))'*diag(risk_weights)*X(:,3);
    
    real_q(3) = 1 - normcdf(p_1,0,sqrt(var_1))- normcdf(p_2,0,sqrt(var_2))+mvncdf([p_1 p_2],[0,0],[var_1 cov_1;cov_1 var_2]+1e-6*eye(2));
    
    % at the end of fourth step 
    t11 = sqrt(tau)*vol_matrix*diag(k)*X(:,4) + sqrt(tau)*vol_matrix*(S(:,2)+S(:,3)+S(:,4));
    t12 = sqrt(tau)*vol_matrix*diag(k)*X(:,4) + sqrt(tau)*vol_matrix*(S(:,3)+S(:,4));
    t13 = sqrt(tau)*vol_matrix*diag(k)*X(:,4)+ sqrt(tau)*vol_matrix*S(:,4);
    t14 = sqrt(tau)*vol_matrix*diag(k)*X(:,4);
    t1 = [t11;t12;t13;t14];
    var_1 = dot(t1,t1);
    q11 = sqrt(tau)*vol_matrix*diag(risk_weights)*X(:,4);
    q1 = [q11;q11;q11;q11];
    var_2 = dot(q1,q1);
    cov_1 = dot(t1,q1);
    %coordinates of the point we want to evaluate
    p_1 = l-(initial_price-beta*S(:,1))'*S(:,1) - (initial_price - gamma*S(:,1)-beta*S(:,2))'*S(:,2)-...
        (initial_price - gamma*(S(:,1)+S(:,2))-beta*S(:,3))'*S(:,3)-...
        (initial_price - gamma*(S(:,1)+S(:,2)+S(:,3))-beta*S(:,4))'*S(:,4)-...
        (initial_price-gamma*(S(:,1) + S(:,2)+S(:,3)+S(:,4)))'*diag(k)*X(:,4);
    p_2 = -(initial_price - gamma*(S(:,1)+S(:,2)+S(:,3)+S(:,4)))'*diag(risk_weights)*X(:,4);
    
    real_q(4) = 1 - normcdf(p_1,0,sqrt(var_1))- normcdf(p_2,0,sqrt(var_2))+mvncdf([p_1 p_2],[0,0],[var_1 cov_1;cov_1 var_2]+1e-6*eye(2));
 
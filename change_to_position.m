function [position] = change_to_position(X_0,strats)
    strats = horzcat(zeros(size(strats,1), 1), strats);
    position = X_0*ones(size(strats,1),1) - cumsum(strats,2);
    

function  plotPortfolioPosition(sol_nonstraint,opt_sol,opt_special,index)
    nonConsPosition = change_to_position(305,sol_nonstraint);
    consPosition = change_to_position(305,opt_sol);
    affine_position = change_to_position(305,opt_special);
    plot([0:5],nonConsPosition(index,:),'-k');
    hold on;
    plot([0:5],consPosition(index,:),'-b');
    plot([0:5],affine_position(index,:),'-r');
    % Set the x-axis to show only integer tick marks
    xticks(0:1:5);
    xticklabels(0:1:5);
    ylabel('remaining position')
    legend('uncosntrained','optimal strategy','optimal affine');
    xlabel('time')
    %title('Liquidation trajactory (risk aversion = 0.1)')
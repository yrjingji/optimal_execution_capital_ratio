
consBigRatioPosition = change_to_position(305,consBigRatio);
plot([0:5],nonConstraintPosition(end,:),'-k');
hold on;
plot([0:5],consSmallRatioPosition(end,:),'-b');
plot([0:5],consBigRatioPosition(end,:),'-m');
% Set the x-axis to show only integer tick marks
xticks(0:1:5);
xticklabels(0:1:5);
ylabel('remaining position')
legend('unconstrained','capital ratio: 6.25% ','capital ratio: 8.75% ');
xlabel('time')
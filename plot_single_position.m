
consBigRatioPosition = change_to_position(305,consBigRatio);
plot([0:5],nonConstraintPosition(11,:),'-k');
hold on;
plot([0:5],consSmallRatioPosition(11,:),'-b');
plot([0:5],consBigRatioPosition(11,:),'-r');
% Set the x-axis to show only integer tick marks
xticks(0:1:5);
xticklabels(0:1:5);
ylabel('remaining position')
legend('optimal trading strategy without chance constraint','capital ratio requirement 6.25% ','capital ratio requirement 8.75% ');
xlabel('time point')
title('single asset: max mean-0.1var')
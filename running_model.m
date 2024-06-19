% running Danino model without coupling for 500 minutes
danino(0, 0, 0, [0, 500], 1)
%% 
% assigning values
AiiA = ans.y(1, :); %#ok<*NOANS>
LuxI = ans.y(2, :);
AHLi = ans.y(3, :);
AHLe = ans.y(4, :);
time = ans.x;

%%
% plotting
% figure(1);
% hold on
% plot(time, AiiA)
% plot(time, LuxI)
% legend("AiiA", "LuxI")
% hold off
% 
% figure(2);
% hold on
% plot(time, AHLi)
% plot(time, AHLe)
% legend("Internal AHL", "External AHL")
% hold off

%% 
% what's the smallest time step?
tmin = min(diff(time));

% an array with evenly spaced time points, based on the smallest time step
timeminn = 0:tmin:max(time);

% for now I will interpolate LuxI and imagine that it is same as GFP which
% I would observe
LuxIminn = interp1(time,LuxI, timeminn, "spline");
% % graph
% figure(4);
% plot(time,LuxI,'o',timeminn,LuxIminn,':.');

%% 
% compute period based on LuxI
LuxIperiod = compute_period(LuxIminn,tmin);



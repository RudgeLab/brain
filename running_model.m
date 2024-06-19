% running Danino model without coupling for 500 minutes
danino(0, 0, 0, [0, 500], 0)
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
% xlabel("Time(min)")
% ylabel("AU")
% hold off
% 
% figure(2);
% hold on
% plot(time, AHLi)
% plot(time, AHLe)
% legend("Internal AHL", "External AHL")
% xlabel("Time(min)")
% ylabel("AU")
% hold off

%% 
tmin = min(diff(time)); % smallest time step
timeminn = 0:tmin:max(time); % an array with evenly spaced time points
LuxIminn = interp1(time,LuxI, timeminn, "spline"); % interpolare LuxI

% graph
% figure(4);
% plot(time,LuxI,'o',timeminn,LuxIminn,':.');
% legend("LuxI", "Intrapolated LuxI")
% xlabel("Time(min)")
% ylabel("AU")

LuxIperiod = compute_period(LuxIminn,tmin); % compute period based on LuxI



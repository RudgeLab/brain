% running Danino model without coupling for 500 minutes
danino(0, 0, 0, [0, 500], 1)
%% 
% assigning values
AiiA = ans.y(1, :);
LuxI = ans.y(2, :);
AHLi = ans.y(3, :);
AHLe = ans.y(4, :);
time = ans.x;

%%
% plotting
figure(1);
hold on
plot(time, AiiA)
plot(time, LuxI)
legend("AiiA", "LuxI")
hold off

figure(2);
hold on
plot(time, AHLi)
plot(time, AHLe)
legend("Internal AHL", "External AHL")
hold off
%%
% plotting
indexoftime = 1:length(time);
figure(3);
hold on
plot(indexoftime, time)
legend("Time")
hold off
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
% plotting time points
indexoftime = 1:length(time);
figure(3);
hold on
plot(indexoftime, time)
legend("Time")
hold off
%% 
% interpolation of the dataset
listofdiff = zeros(2, 1);
for i=1:length(time)
    if i>1
        difference = time(i)-time(i-1); % difference between current and previous time point
        listofdiff = cat(2, listofdiff, [i;difference]); % add which time point it is and how big is the difference
    end
end
listofdiff(2,1) = NaN; % ensure that first position is empty
% what's the smallest time step?
minn = min(listofdiff(2,:));

% an array with evenly spaced time points, based on the smallest time step
timeminn = 0:minn:500;

% for now I will interpolate LuxI and imagine that it is same as GFP which
% I would observe
LuxIminn = interp1(time,LuxI, timeminn);
% graph
plot(time,LuxI,'o',timeminn,LuxIminn,':.');
%% 
% % alternative array with measurement every 10 mins, to replicate measuring
% % system every 10 minutes
% timetenm = 0:10:500;
% % for now I will interpolate LuxI and imagine that it is same as GFP which
% % I would observe
% LuxItenm = interp1(time,LuxI, timetenm);
% % graph
% plot(time,LuxI,':.',timetenm,LuxItenm,'o'); % not a good choice, this
% % misses peaks
%% 



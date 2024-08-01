% running Danino model without coupling for 500 minutes
danino(0, 0, 0, [0, 500], 0)

% assigning values
AiiA = ans.y(1, :); %#ok<*NOANS>
LuxI = ans.y(2, :);
AHLi = ans.y(3, :);
AHLe = ans.y(4, :);
time = ans.x;

t = 1; % time step
timeq = 0:t:max(time); % an array with evenly spaced time points equal to t
LuxIeq = interp1(time,LuxI, timeq, "spline"); % interpolare LuxI
LuxIperiod = compute_period(LuxIeq,t); % compute period based on LuxI 

%% 
couplings = logspace(-9, 0, 10);
% couplings = linspace(0.0001, 0.001, 10);
SHIL_coupling(@danino, [0, 500], couplings, 2)
%% 

% SHIL(@danino, [0, 500], 0.0001, 2)
%% 
% code to save figures - PC specific
FolderName = "C:\Users\Luiza\Desktop\Uni_code\Danino exp\SHIL_coupling";   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  saveas(FigHandle,fullfile(FolderName, [FigName '.png'])); %Specify format for the figure
end

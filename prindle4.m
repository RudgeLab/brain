ntrain = 5000;
load ('mnist.mat');
XTrain = training.images(:,:,1:ntrain);
YTrain = training.labels(1:ntrain) + 1;
%[XTrain,YTrain,anglesTrain] = digitTrain4DArrayData;
%[XTest,YTest,anglesTest] = digitTest4DArrayData;

nx = 32;
ny = 32;
nreps = 50;
for i=1:nreps
	% Use node id as random seed for selecting training image
	%stream = RandStream('mt19937ar','seed', (input_idx-1) * 5 + i - 1);
	%RandStream.setGlobalStream(stream);

	ntrain = size(XTrain, 3);
	zero_idxs = find(double(YTrain)==1);
	one_idxs = find(double(YTrain)==2);
	all_idxs = [zero_idxs; one_idxs];
	%indices = ceil(rand(1,1)*numel(all_idxs));
	%idx = all_idxs((input_idx - 1) * nreps + i);
	idx = (input_idx - 1) * nreps + i;

	% Use same seed for all nodes to generate same mapping
	stream = RandStream('mt19937ar','seed', 1);
	RandStream.setGlobalStream(stream);

	mapping = randperm(nx*ny, 28*28);
	input_set = zeros(nx,ny);
	input_phase = zeros(nx,ny);
	input_set(mapping) = 1;
	input_phase(mapping) = XTrain(:,:,idx);

	tic
	[time_pts, im] = run_sim(input_phase, input_set);
	toc

	t = time_pts;
	tmin = max(t)*0.75;
	tmax = max(t)*1;
	tt = tmin:1:tmax;
	
	%sim = zeros(nx,nx,length(tt));
	%for x=1:nx
	%	for y=1:ny
	%		y1 = interp1(time_pts, reshape(im(x,y,1,:), [], 1), tt);		
	%		sim(x,y,:) = y1;
	%	end
	%end
	%save(sprintf("sim_%d.mat", idx), "sim");

	[pks,locs] = findpeaks(reshape(im(1,1,1,:), [], 1));
	period= mean(diff(t(locs)));

	phases = zeros(nx,ny);
	%y0 = interp1(t, reshape(im(1,1,1,:), [], 1), tt);
	T = 50.1781;
	y0 = 0.5 * (1 + sin(2 * pi * tt / T));
	[nx,ny,nc,nt] = size(im);
	for x=1:nx
	    for y=1:ny
		y1 = interp1(t, reshape(im(x,y,1,:), [], 1), tt);
		[c,lag]=xcorr(y0, y1);
		[maxC,I]=max(c);
		lag = lag(I);
		p = mod(lag * 360 / period, 360);
		if p>180
		    p = 360 - p;
		end
		% phases=1 in phase, phase=-1 antiphase
		phases(x,y) = 1 - 2 * (p / 180);
	    end
	end

	save(sprintf("phases_%d.mat", idx), "phases");
end

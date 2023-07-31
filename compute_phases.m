function phases = compute_phases(time_pts, im)
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
	T = 50.1781; % period
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

end
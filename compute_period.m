function period = compute_period(s, dt) % s is the signal, dt is the length of the time step
    ac = xcorr(s, s);
    [~,locs]=findpeaks(ac); % find indexes of the peaks in ac
    period = mean(diff(locs(2:end-1))) * dt; % calculate the index differences between peaks (with the exlusion  of 1st and last), find mean and multiply by the time step
end
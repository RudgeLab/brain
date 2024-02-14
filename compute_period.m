function period = compute_period(s, dt)
    ac = xcorr(s, s);
    [~,locs]=findpeaks(ac);
    period = mean(diff(locs(2:end-1))) * dt;
end
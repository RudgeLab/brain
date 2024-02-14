% compute phase between 2 signals
function p = compute_phase(s1, s2, period, dt)
    [cc,lags] = xcorr(s1, s2);
    ic = round(length(cc)/2);
    cc = cc(ic:end);
    lags = lags(ic:end);
    [pks,locs] = findpeaks(cc);
    lag = lags(locs(1));
    p = 360 * lag / period * dt;
    % plot(360 * lags * dt / period, cc)
    if p>180
        p = 360 - p;
    end
end

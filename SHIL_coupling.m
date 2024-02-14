% Demonstrate injection locking of Danino oscillator to input AHL
function periods = SHIL_coupling(model, tspan, couplings, n)
    ts = round(mean(tspan));
    T = linspace(tspan(1), tspan(2), 1440);
    dt = mean(diff(T));
    phase = 0
    periods = []

    % Reference with no coupling to compute period
    sol = model(1, 0, 0, tspan, n);
    y = sol.y;
    t = sol.x;
    iy1ref = interp1(t, y(1,:), T);
    plot(T, iy1ref);
    % Compute period from last half of data
    period  = compute_period(iy1ref, dt);
    period
    
    for coupling = couplings
        input_signal = prof_pulse(T, period/n, phase*n);
        ref_signal = prof_pulse(T, period, phase);
        
        % Simulate system with input at SHIL frequency (half natural period)
        sol = model(period, phase, coupling, tspan, n);
        y = sol.y;
        t = sol.x;
        iy1 = interp1(t, y(1,:), T);

        p = compute_period(iy1, dt)
        periods(end+1) = p;
    
        figure();
        hold on;
        plot(T(ts:end), ref_signal(ts:end) * max(iy1), 'r--');
        plot(T(ts:end), input_signal(ts:end) * max(iy1), 'g--');
        plot(T(ts:end), iy1(ts:end), 'b');
        title(gca, sprintf('%f', coupling));
    end
    
    figure();
    plot(couplings,  periods,  'r.', markersize=20);
    xlabel('coupling');
    ylabel('period');
end
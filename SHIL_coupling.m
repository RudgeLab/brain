% Demonstrate injection locking of Danino oscillator to input AHL
function periods = SHIL_coupling(model, tspan, couplings, n)
% model - input @_model_
% tspan - time range in mins
% couplings - a range of how much of SHIL molecule is added per time step
% n - is how many times higher is the frequency of SHIL signal compared to
                                                   % the natural frequency
    ts = round(mean(tspan));
    T = linspace(tspan(1), tspan(2), 1440); % generate 1440 evenly spaced time points
    dt = mean(diff(T));
    phase = 0
    periods = []

    % Reference with no coupling to compute period
    sol = model(1, 0, 0, tspan, n);
    y = sol.y;
    t = sol.x;
    iy1ref = interp1(t, y(1,:), T); % interpolate AiiA, non-spline
    plot(T, iy1ref); % plot AiiA against time
    legend("AiiA")
    xlabel("Time(t)")
    ylabel("AU")
    title("Without coupling")
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
        iy1 = interp1(t, y(1,:), T); % interpolate AiiA, non-spline

        p = compute_period(iy1, dt)
        periods(end+1) = p;
    
        figure();
        hold on;
        plot(T(ts:end), ref_signal(ts:end) * max(iy1), 'r--');
        plot(T(ts:end), input_signal(ts:end) * max(iy1), 'g--');
        plot(T(ts:end), iy1(ts:end), 'b');
        title(sprintf('Coupling = %0.3g', coupling));
        legend("?", "SHIL signal", "AiiA")
        xlabel("Time(t)")
        ylabel("AU")

    end
    
    figure();
    if std(diff(couplings))>=0.000001 % check standard deviation of couplings 
                                 % to determine if x is log or linear
        semilogx(couplings,  periods,  'r.', markersize=20)
    else 
        plot(couplings,  periods,  'r.', markersize=20);
    end
    xlabel('Coupling (AU)');
    ylabel('Period(t)');
    title("Correlation between couplings and period")
    grid on
end
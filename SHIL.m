% Demonstrate injection locking of Danino oscillator to input AHL
function [phases, ref_phase_lags, phase_lags, shil_phase_lags] = SHIL(model, tspan, coupling, n)
    T = linspace(tspan(1), tspan(2), 1440);
    dt = mean(diff(T));
    
    % Reference with no coupling to compute period
    sol = model(68, 0, 0, tspan, n);
    y = sol.y;
    t = sol.x;
    iy1ref = interp1(t, y(1,:), T);
    plot(T, iy1ref);
    % Compute period from last half of data
    period  = compute_period(iy1ref(721:end), dt);
    period

    ts = round(1440 - period * 6);  
    
    phase_lags = [];
    ref_phase_lags = [];
    shil_phase_lags = [];
    % Input and reference signal
    phases = linspace(0, pi, 8);
    for phase = phases
        input_signal = prof_pulse(T, period/n, phase*n);
        ref_signal = prof_pulse(T, period, phase);
        
        % Simulate system with input at SHIL frequency (half natural period)
        sol = model(period, phase, coupling, tspan, n);
        y = sol.y;
        t = sol.x;
        iy1 = interp1(t, y(1,:), T);
        figure();
        plot(T, iy1);

        % Compute phase lag of result to reference signal, using last half of data
        phase_lag = compute_phase(iy1(ts:end), ref_signal(ts:end), period, dt)
        phase_lags(end+1) = phase_lag;
    
        % Compute phase lag of reference oscillator to ref signal, using last half of data
        ref_phase_lag = compute_phase(iy1ref(ts:end), ref_signal(ts:end), period, dt)
        ref_phase_lags(end+1) = ref_phase_lag;

        % Compute phase lag between SHIL and no SHIL
        shil_phase_lag = compute_phase(iy1ref(ts:end), iy1(ts:end), period, dt)
        shil_phase_lags(end+1) = shil_phase_lag;
        %figure();
        %plot(lags, cc)
    
        figure();
        hold on;
        plot(T(ts:end), ref_signal(ts:end) * max(iy1), 'r--');
        plot(T(ts:end), input_signal(ts:end) * max(iy1), 'g--');
        plot(T(ts:end), iy1(ts:end), 'b');
        plot(T(ts:end), iy1ref(ts:end), 'k');
        title(gca, sprintf('%f', phase_lag));
    end
    
    ref_phase_lags
    phase_lags
    figure();
    plot(ref_phase_lags,  phase_lags,  'r.', markersize=20);
    hold on;
    plot(phases*180/pi,  phase_lags,  'g.', markersize=20);
    plot([90,90], [0,180], 'k--');
    plot([100,100], [0,180], 'k--');
    xlabel('Initial phase (degs)');
    ylabel('Final phase (degs)');
end
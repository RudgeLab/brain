
function sol = danino(period, phase, coupling, tspan, n)
    lags = [10];        % time delays
    CA = 1;             % Copy number of AiiA
    CI = 4;             % Copy number of LuxI
    del = 1e-3;         % Promoter leakiness
    alpha = 2500;       % Promoter strength
    tau = 10;           % time delay for the delay differential equation
    k = 1;              % derived from switching concentration
    k1 = 0.1;           % derived from switching concentration
    b = 0.06;           % AHL synthesis rate by LuxI
    gammaA = 15;        % degradation of AiiA via proteases rate
    gammaI = 24;        % degradation of LuxI via proteases rate
    gammaH = 0.01;      % degradation of AHL via AiiA rate
    f = 0.3;            % derived from switching concentration
    g = 0.01;           % derived from switching concentration
    d0 = 0.88;          % related to slowing down due to the cell density
    D = 2.5;            % Diffusion rate/constant
    mu = 0.6;           % Dilution via flow
    d = 0.7;            % cell density
       
    sol = solve(lags, tspan, CA, CI, del, alpha, k, k1, b, gammaA, gammaI, gammaH, f, g, d, d0, D, mu, period, phase, coupling, n);
    y = sol.y;
    time_pts = sol.x;
end

function sol = solve(lags, tspan, CA, CI, del, alpha, k, k1, b, gammaA, gammaI, gammaH, f, g, d, d0, D, mu, period, phase, coupling, n)
    sol = dde23(@ddefun, lags, @history, tspan); 
    
    function dydt = ddefun(t,y,Z)
      Hlag = Z(3,1);    % Past concentration of Internal AHL
      A = y(1);         % AiiA
      I = y(2);         % LuxI
      Hi = y(3);        % Internal AHL
      He = y(4);        % External AHL

      % Addition of external AHL
      Hetot = He + coupling * shil_signal(t, n);                        % Total external AHL after addition of SHIL signal

      P = (del + alpha*Hlag^2) / (1 + k1*Hlag^2);                       % Hill function, protein expression

      dAdt = CA * (1 - (d/d0)^4) * P - gammaA * A / (1 + f*(A+I));      % change in AiiA
      dIdt = CI * (1 - (d/d0)^4) * P - gammaI * I / (1 + f*(A+I));      % change in LuxI
      dHidt = b*I/(1 + k*I) - gammaH*A*Hi / (1 + g*A) + D*(Hetot-Hi);   % change in internal AHL
      dHedt = -d / (1 - d) * D*(He-Hi) - mu*He;                         % change in external AHL

      dydt = [dAdt; dIdt; dHidt; dHedt];                                % compile into an array
    end

    function s = history(t)
        s = zeros(4,1);                                                 % make an array of zeros
        s(2) = 100;                                                     % Set initial LuxI to a 100 units
    end

    function s = shil_signal(t, n)
        s = prof_pulse(t, period/n, phase*n); 
        % if t<period*4
        %     s = prof_pulse(t, period, phase);
        % elseif t<period*8
        %     s = 0.5 * (prof_pulse(t, period, phase) + prof_pulse(t, period/n, phase*n));
        % else
        %     s = prof_pulse(t, period/n, phase*n);
        % end
    end
end





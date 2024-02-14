
function sol = danino(period, phase, coupling, tspan, n)
    lags = [10];
    CA = 1;
    CI = 4;
    del = 1e-3;
    alpha = 2500;
    tau = 10;
    k = 1;
    k1 = 0.1;
    b = 0.06;
    gammaA = 15;
    gammaI = 24;
    gammaH = 0.01;
    f = 0.3;
    g = 0.01;
    d0 = 0.88;
    D = 2.5;
    mu = 0.6;  
    d = 0.7;
       
    sol = solve(lags, tspan, CA, CI, del, alpha, k, k1, b, gammaA, gammaI, gammaH, f, g, d, d0, D, mu, period, phase, coupling, n);
    y = sol.y;
    time_pts = sol.x;
end

function sol = solve(lags, tspan, CA, CI, del, alpha, k, k1, b, gammaA, gammaI, gammaH, f, g, d, d0, D, mu, period, phase, coupling, n)
    sol = dde23(@ddefun, lags, @history, tspan);
    
    function dydt = ddefun(t,y,Z)
      Hlag = Z(3,1);
      A = y(1);
      I = y(2);
      Hi = y(3);
      He = y(4);

      % Addition of external AHL
      Hetot = He + coupling * shil_signal(t, n);

      P = (del + alpha*Hlag^2) / (1 + k1*Hlag^2);

      dAdt = CA * (1 - (d/d0)^4) * P - gammaA * A / (1 + f*(A+I));
      dIdt = CI * (1 - (d/d0)^4) * P - gammaI * I / (1 + f*(A+I));
      dHidt = b*I/(1 + k*I) - gammaH*A*Hi / (1 + g*A) + D*(Hetot-Hi);
      dHedt = -d / (1 - d) * D*(He-Hi) - mu*He;

      dydt = [dAdt; dIdt; dHidt; dHedt];
    end

    function s = history(t)
        s = zeros(4,1);
        s(2) = 100;
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





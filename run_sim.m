
function [time_pts, im] = run_sim(input_phase, input_set)
    lags = [10 20];
    nx = 32;
    ny = 32;
    alpha = 8.25;
    v = 1;
    C0 = 6;
    gamma = 5.75;
    k = 10/4;
    mu = 20;
    alpha_p = 1;
    gamma_p = 10;
    D = 7;
    delta = 2;
    
    tspan = [0 1000];
    
    sol = solve(lags, tspan, nx, ny, alpha, v, C0, gamma, k, mu, alpha_p, gamma_p, D, delta, input_phase, input_set);
    im = reshape(sol.y, nx, ny, 2, []);
    time_pts = sol.x;
end

function sol = solve(lags, tspan, nx, ny, alpha, v, C0, gamma, k, mu, alpha_p, gamma_p, D, delta, input_phase, input_set)
    mag = 50;
    T = 50.1781;
    sol = dde23(@ddefun, lags, @history, tspan);
    
    function dydt = ddefun(t,y,Z)
      ylag1 = reshape(Z(:,1), nx, ny, 2);
      Xlag = ylag1(:,:,1);
      ylag2 = reshape(Z(:,2), nx, ny, 2);
      Plag = ylag2(:,:,2);
      y = reshape(y, nx, ny, 2);
      X = y(:,:,1);
      P = y(:,:,2);

      dXdt = zeros(size(X));
      dPdt = zeros(size(P));
      for x=1:nx;
          for y=1:ny;
              delP = - 4 * P(x,y);
              if x>1;
                  delP = delP + P(x-1,y);
              else
                  delP = delP + P(end,y);
              end
              if x<nx;
                  delP = delP + P(x+1,y);
              else
                  delP = delP + P(1,y);
              end
              if y>1;
                  delP = delP + P(x,y-1);
              else
                  delP = delP + P(x,end);
              end
              if y<ny;
                  delP = delP + P(x,y+1);
              else
                  delP = delP + P(x,1);
              end
              dXdt(x,y) = alpha * (1 + v*Plag(x,y)) / ((1 + Xlag(x,y)/C0)^2) - gamma * X(x,y) / (k + X(x,y));
              dPdt(x,y) = mu + alpha_p * X(x,y) - gamma_p * P(x,y) + D * delP / (delta^2);
              if input_set(x,y)
                  dPdt(x,y) = dPdt(x,y) + input(t, mag, T, input_phase(x,y)*pi/180);
              end
          end
      end

      dydt = [reshape(dXdt, nx*ny, 1); reshape(dPdt, nx*ny, 1)];
    end

    function s = history(t)
      s = zeros(nx*ny*2, 1);
    end

    function y = input(t, mag, T, phi)
      y = mag * 0.5 * (1 + sin(2 * pi * t / T + phi));
    end
    
    % function y = input(t, mag, T, phi)
    %     y = mag * (mod(t + phi*T/2/pi, T) > T/2);
    % end

end





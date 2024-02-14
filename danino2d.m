
function sol = danino2d(nx, ny, dx, tspan)
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
    D1 = 2000;
       
    sol = solve(lags, tspan, CA, CI, del, alpha, k, k1, b, gammaA, gammaI, gammaH, f, g, d, d0, D, mu, D1, nx, ny, dx);
    y = sol.y;
    time_pts = sol.x;
end

function sol = solve(lags, tspan, CA, CI, del, alpha, k, k1, b, gammaA, gammaI, gammaH, f, g, d, d0, D, mu, D1, nx, ny, dx)
    sol = dde23(@ddefun, lags, @history, tspan);
    
    function dydt = ddefun(t,y,Z)
      y = reshape(y,[4,nx,ny]);
      dydt = zeros(4,nx,ny);
      ylag = reshape(Z,[4,nx,ny]);
      Hlag = reshape(ylag(3,:,:), nx,ny);
      A = reshape(y(1,:,:), nx,ny);
      I = reshape(y(2,:,:), nx,ny);
      Hi = reshape(y(3,:,:), nx,ny);
      He = reshape(y(4,:,:), nx,ny);

      for i=[1:nx]
          for j=[1:ny]
              delHe = - 4 * He(i,j);
              if i>1;
                  delHe = delHe + He(i-1,j);
              else
                  delHe = delHe + He(end,j);
              end
              if i<nx;
                  delHe = delHe + He(i+1,j);
              else
                  delHe = delHe + He(1,j);
              end
              if j>1;
                  delHe = delHe + He(i,j-1);
              else
                  delHe = delHe + He(i,end);
              end
              if j<ny;
                  delHe = delHe + He(i,j+1);
              else
                  delHe = delHe + He(i,1);
              end
              delHe = delHe / 4;

              P = (del + alpha*Hlag(i,j)^2) / (1 + k1*Hlag(i,j)^2);
        
              dAdt = CA * (1 - (d/d0)^4) * P - gammaA * A(i,j) / (1 + f*(A(i,j)+I(i,j)));
              dIdt = CI * (1 - (d/d0)^4) * P - gammaI * I(i,j) / (1 + f*(A(i,j)+I(i,j)));
              dHidt = b*I(i,j)/(1 + k*I(i,j)) - gammaH*A(i,j)*Hi(i,j) / (1 + g*A(i,j)) + D*(He(i,j)-Hi(i,j));
              dHedt = -d / (1 - d) * D*(He(i,j)-Hi(i,j)) - mu*He(i,j) + D1 * delHe / (dx^2);
        
              dydt(:,i,j) = [dAdt; dIdt; dHidt; dHedt];
          end
      end
      dydt = reshape(dydt, [], 1);
    end
    function s = history(t)
        s = zeros(4,nx,ny);
        s(2,:,:) = 100  + 100 * rand(nx,ny);
    end
end





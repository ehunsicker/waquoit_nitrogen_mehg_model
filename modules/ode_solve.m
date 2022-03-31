function [ts,yy] = ode_solve(odefun, ts, y0)
% ODE_SOLVE
%   Solve ODEs using Euler method to handle stochastic parameters 
yy = zeros(numel(ts),numel(y0));
yy(1,:) = y0;
for i=2:numel(ts)
    clear dydt
    dt = ts(i)-ts(i-1);
    dydt = odefun(ts(i-1),yy(i-1,:));
    yy(i,:) = yy(i-1,:) + dydt.'*dt;
    yy(i,:) = max(yy(i,:),y0);
end
end


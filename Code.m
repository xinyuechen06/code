clc;
clear;
N = 5; % Number of nodes
dt = 0.01; % Time step
steps = 300000; % Integration steps
omega = 1; % Natural frequency
a = 0.5;
alpha = 0.5;
epsilon1 = 1;
epsilon2 = 2;
epsilon = [epsilon1; epsilon2*ones(N-1,1)];

A = [0 ones(1,N-1); ones(N-1,1) zeros(N-1)];
G = alpha*diag(sum(A, 2)) - A;
G(1,:) = G(1,:) / (N-1);

% Initial conditions
xold = 2*rand(N, 1) - 1;
yold = 2*rand(N, 1) - 1;
x = zeros(steps, N);
y = zeros(steps, N);

% Numerical integration
for t = 1:steps
    couplingx = -G*xold;
    couplingy = -G*yold;
    dxdt = a*xold - omega*yold - xold.*(xold.^2 + yold.^2) + epsilon.*couplingx;
    dydt = omega*xold + a*yold - yold.*(xold.^2 + yold.^2) + epsilon.*couplingy;
    xnew = xold + dt*dxdt;
    ynew = yold + dt*dydt;
    xold = xnew;
    yold = ynew;
    x(t, :) = xnew;
    y(t, :) = ynew;
end

x = round(x, 4);
y = round(y, 4);
theta = atan2(y, x);
rho = sqrt(y.^2 + x.^2);
theta(theta < 0) = theta(theta < 0) + 2*pi;

% Calculate new metrics tau
thetamatrix = zeros(N, N);
for ii = 1:N
    for jj = ii+1:N
        thetamatrix(ii, jj) = theta(end, jj) - theta(end, ii);
    end
end
cosmatrix = triu(cos(thetamatrix), 1);
tau_direct = round(mean(cosmatrix(1, 2:N)), 4);
tau_indirect = round(2*sum(sum(cosmatrix(2:N, 2:N))) / ((N-1) * (N-2)), 4);

% Calculate old metrics r
errormatrix = zeros(N, N);
for ii = 1:N
    for jj = ii+1:N
        errormatrix(ii, jj) = abs(mean(exp(1i * (theta(end/2:end, ii) - theta(end/2:end, jj)))));
    end
end
r_direct = mean(errormatrix(1, 2:N));
r_indirect = sum(sum(triu(errormatrix(2:end, 2:end), 1))) * 2 / (N-1) / (N-2);



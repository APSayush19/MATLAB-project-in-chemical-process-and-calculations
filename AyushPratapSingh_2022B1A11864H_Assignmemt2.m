% Ayush_2022b1A11864H_Assignment2.m

clc;
clear;

% Given Parameters
h = 100; % W/m^2K
k = 180; % W/mK
d = 0.005; % m
L = 0.05; % m
T_inf = 25; % deg C
T_b = 100; % deg C

% Derived Parameters
Ac = pi * d^2 / 4;
P = pi * d;
Af = pi * d * L;
m = sqrt(h * P / (k * Ac));

% User input: number of intervals
n = input('Enter number of intervals (n): ');
dx = L / n;
x = linspace(0, L, n+1)'; % column vector

% Coefficient for numerical scheme
a = h * P / (k * Ac);

% Initialize matrix A and vector b
A = zeros(n+1);
b = zeros(n+1, 1);

% Boundary condition at x = 0 (base): T(0) = T_b
A(1, 1) = 1;
b(1) = T_b;

% Interior nodes: central difference
for i = 2:n
    A(i, i-1) = 1/dx^2;
    A(i, i) = -2/dx^2 - a;
    A(i, i+1) = 1/dx^2;
    b(i) = -a * T_inf;
end

% Tip condition at x = L (convective tip)
A(n+1, n) = -k / dx;
A(n+1, n+1) = k / dx + h;
b(n+1) = h * T_inf;

% Solve Ax = b using Gauss Elimination
T_numeric = gauss_elimination(A, b);

% Analytical solution
theta_b = T_b - T_inf;
theta = @(x) (cosh(m*(L - x)) + (h/(m*k)) * sinh(m*(L - x))) ./ ...
             (cosh(m*L) + (h/(m*k)) * sinh(m*L));
T_analytical = T_inf + theta(x) * theta_b;

% Plotting
figure;
plot(x, T_numeric, 'ro-', 'LineWidth', 2);
hold on;
plot(x, T_analytical, 'b--', 'LineWidth', 2);
xlabel('Position along fin, x (m)');
ylabel('Temperature (°C)');
title('Temperature Distribution along Fin');
legend('Numerical Solution', 'Analytical Solution');
grid on;

% Heat transfer rate at base
q_f = -k * Ac * (T_numeric(2) - T_numeric(1)) / dx;

% Fin effectiveness
effectiveness = q_f / (h * Ac * (T_b - T_inf));

% Fin efficiency
efficiency = q_f / (h * Af * (T_b - T_inf));

% Display results
fprintf('\nTemperature at tip of the fin: %.4f °C\n', T_numeric(end));
fprintf('Fin Effectiveness: %.4f\n', effectiveness);
fprintf('Fin Efficiency: %.4f\n', efficiency);


% ------------------------
% GAUSS ELIMINATION METHOD
% ------------------------
function x = gauss_elimination(A, b)
    n = length(b);
    Aug = [A b];

    % Forward Elimination
    for i = 1:n-1
        for j = i+1:n
            factor = Aug(j,i)/Aug(i,i);
            Aug(j,:) = Aug(j,:) - factor * Aug(i,:);
        end
    end

    % Back Substitution
    x = zeros(n,1);
    x(n) = Aug(n,end)/Aug(n,n);
    for i = n-1:-1:1
        x(i) = (Aug(i,end) - Aug(i,i+1:n) * x(i+1:n)) / Aug(i,i);
    end
end

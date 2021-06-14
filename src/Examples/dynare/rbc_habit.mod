% Neoclassical Investment Model, Log Utility with External Habit Formation

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y c k z s;
varexo e;


%----------------------------------------------------------------
% 2. Parameters values
%----------------------------------------------------------------

parameters beta delta alpha rho phi sigma;
alpha   = 0.6;
beta    = 0.96;
delta   = 0.10;
rho     = 0.95;
phi     = 0.05;
sigma   = 0.01;

%----------------------------------------------------------------
% 3. Model (Euler Equation, Budget Constraint, Laws of Motion)
%----------------------------------------------------------------

model;
  (c-phi*s)^(-1) = beta*(c(+1)-phi*s(+1))^(-1)*(alpha*exp(z(+1))*(k^(alpha-1))+1-delta); % Euler
  c+k-(1-delta)*k(-1) = y; % Budget Constraint
  s = c(-1); % Habit Stock
  y = exp(z)*(k(-1)^alpha); % Production Function
  z = rho*z(-1)+e; % TFP AR(1)
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

%%Steady State (Educated) Guesses

initval;

  k = ((1/beta-1+delta)/alpha)^(1/(alpha-1));
  c = k^alpha-k;
  s = c;
  z = 0;
  e = 0;

end;

shocks;
var e = sigma^2;
end;

steady;

stoch_simul;

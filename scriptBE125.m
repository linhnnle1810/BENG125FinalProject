%%%%%%% BENG 125 Final Project 2019 %%%%%%%%
%%%%%%% Team Members %%%%%%%
% Jose Zaragoza, Allen Le, Jacqueline Olness, Kokila Perera, and Linh N Le %
% Date: 5/2019 %

%% Description of project:
% To model the trends of meme popularity using models defined in 
%Wang, L., & Wood, B. C. (2011). An epidemiological approach to model the viral propagation of memes. 
%Applied Mathematical Modelling, 35(11), 5442-5447. doi:10.1016/j.apm.2011.04.035

% Pseudo-code
% Preparing the model:
% 1. Define differential equations with unknown/unfitted parameters
%    Susceptible population: dS/dt = -alpha*S*I
%    Infected population: dI/dt = (alpha - beta)*S*I + (beta*N - gamma -
%    beta*I)*I
%    Parameters to tune: alpha = transmission rate prop cnst, beta = reinfection
%    rate prop cnst, gamma = loss of interest prop cnst
%    Constants: N = total population
N = 3*10^1;
a = 5.78*10^-1; 
b = 3.91*10^-4;
g = 1.26*10^-2;

[s,i] = meshgrid(-2:0.2:2,-2:0.2:2);
ds = -a*s.*i;
di = (a-b)*s.*i+(b*N-g-b*i)*i;

figure
quiver(s,i,ds,di)
% 2. Linearize the system using Jacobian (as dI/dt is non-linear in terms of I while system is in 2D)
% 3. Find eigenvalues and eigenvectors in terms of tuning parameters, I* &
%    S*
% 4. Load experimental data specific to a meme:
%        - population of infected individuals over time
%
% ***LOOP***
% 5. Define fixed points for the system: E = (S*, I*)
%       For a meme dying out: FP is (Sbar,0) where Sbar is a func of S0 &
%       I0
%       For a meme persisting: FP is (0,N - gamma/beta)
% 6. Solve eigenvalues & vectors for the equilibrium condition, E
% 7. Form equations for S and I as a function of time (using eigenvalues & eigenvectors)
% 8. Use initial conditions to solve for constants 
%    (this will be the final model equation where unknown parameters have to be tuned)
%
% Running modelling: 
% 9. Carry out least squares fitting on model equation with the
%    experimental data
% ***END LOOP***
%
% 10. Plot the infection curve with final tuned parameters
%    (eg. I vs t, or we can modify I as needed to represent some "search term index")
%%%%%%%%%%%%%%%%%%%%%%

%%
clear all; 
clc;
close all; 
format long; 

%% Load experimental data


%% Defining Model

% Define initial tuning parameters
% N = 3e1; time = [0:1:1200]; a0 = 5.78e-1; b0 = 3.91e-4; g0 = 1.26e-2; % parameters are based on 'showdown' values
% N = 1e3; time = [0:1:1400]; a0 = 3.39e-3; b0 = 3.35e-3; g0 = 3.35; % parameters are based on 'O RLY' values
N = 2e2; time = [0:10:20000]; a0 = 1.62e-4; b0 = 1.52e-4; g0 = 3e-2; % parameters are based on 'blog' values

% Define initial conditions
syms alpha beta gamma C0 C1 t S I
t0 = 0; S0 = N - 1; I0 = N - S0;

% Define fixed point
I_f = N - gamma/beta; S_f = 0; % for a persisting meme
% I_f = 0; S_f = 100; % S_f = (gamma-beta*N)/(alpha-beta); % for a dying meme

% Defining Jacobian matrix at FP to obtain eigenvalues and eigenvectors from
J_E = [-1*alpha*I_f, -1*alpha*S_f;...
       (alpha - beta)*I_f, (alpha-beta)*S_f + beta*(N - gamma/beta - 2*I_f)];
   
% Get eigenvalues and eigenvectors
[vec,val] = eig(J_E);
lambda1 = val(1,1); lambda2 = val(2,2);

% Forming time-dependent equations: eqns = [S,I]
eqns = [S == C0*exp(lambda1*t)*vec(1,1) + C1*exp(lambda2*t)*vec(1,2),...
        I == C0*exp(lambda1*t)*vec(2,1) + C1*exp(lambda2*t)*vec(2,2)];

%% Run Model Fitting till least squares fit produces the best result

% Begin LOOP

% Applying tuning parameters
eqns = vpa(subs(eqns,[alpha,beta,gamma],[a0,b0,g0]));

% Apply initial conditions
eqns0 = subs(eqns,[S,I,t],[S0,I0,t0]);

%Solving for constants using initial conditions
[sol_C0, sol_C1] = solve(eqns0,[C0,C1]);
eqns = subs(eqns,[C0,C1],[sol_C0, sol_C1]);

% Obtain values for linear regression
tempS = subs(eqns(1),t,time);
tempI = subs(eqns(2),t,time);

for c = 1:length(time)
    S_final(c) = double(solve(tempS(c),S));
    I_final(c) = double(solve(tempI(c),I));
end

% Carry out least squares fit to change a0, b0, g0

% END LOOP when fitting is sufficient

%% Plot Results

figure(1);
plot(time,I_final);
title('Viral "blog"')
xlabel('days');
ylabel('Infected Population, I(t)');

figure(2);
plot(time,S_final);
title('Viral "blog"')
xlabel('days');
ylabel('Susceptible Population, S(t)');


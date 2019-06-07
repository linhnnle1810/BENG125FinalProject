%%%%%%%%%%%%Parameter Finding Algorithm%%%%%%%%%%%%%%%%%%
clc;
close all;
clear;

%% Initial setup
global views
global time

%Constant Population
N=200;

%Initial Conditions
I0=1;
S0=N-I0;
R0 = 0;

%% Load Data into the Code (Change the file name per example)
A = importdata('GangnamStyle.csv');
Views=A.data(:,2);
Month=A.data(:,1);

%Data imported from google trends: Must be inputted as column vectors
time=Month;
views=Views;

%Initialize parameters with possible values
B=rand(3,1)*0.1;
% B=[5.78e-1,3.91e-4,1.26e-1];

param = lsqnonlin(@g,B);

%% older curve

f=@(t,y) [-B(1)*y(1)*y(2) ; B(1)*y(1)*y(2)- B(2)*y(2)*(y(2)+y(3)); B(2)*y(2)*(y(2)+y(3))];
trange=[0:1:length(time)-1];
initials=[S0 I0 R0];
[t1, infect]=ode45(f,trange,initials);
old_I = infect(:,2);

%% newer curve
B1 = param;

f1=@(t,y) [-B1(1)*y(1)*y(2) ; B1(1)*y(1)*y(2)- B1(2)*y(2)*(y(2)+y(3)); B1(2)*y(2)*(y(2)+y(3))];
trange=[0:1:length(time)-1];
initials=[S0 I0 R0];
[t1, infect]=ode45(f1,trange,initials);
new_I = infect(:,2);

%% Plotting

figure(1);
plot(time,old_I);
hold on;
plot(time,views,'.r');
plot(time,new_I);
title('Gangnam Style')
xlabel('time(months)');
ylabel('number of views');
legend('old fit','data','new fit');
grid on; grid minor;

%% Function definitions

function I = g(B)
    global views
    global time
    
    %Constant Population
    N=200;
    
    %Initial Conditions
    I0=1;
    S0=N-I0;
    R0 = 0;
    
    %Ode Solver
    f=@(t,y) [-B(1)*y(1)*y(2) ; B(1)*y(1)*y(2)- B(2)*y(2)*(y(2)+y(3)); B(2)*y(2)*(y(2)+y(3))];
    trange=[0:1:length(time)-1];
    initials=[S0 I0 R0];
    [t1, infect]=ode45(f,trange,initials);
    I = infect(:,2) - views;
end





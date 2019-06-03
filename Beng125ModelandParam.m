%%%%%%%%%%%%Parameter Finding Algorithm%%%%%%%%%%%%%%%%%%

clc;
close all;
clear;

%Load Data into the Code (Change the file name per example)
A = importdata('GangnamStyle.csv');
Views=A.data(:,2);
Month=A.data(:,1);

%Constant Population
N=200;

%Initial Conditions
I0=1;
S0=N-I0;

%Data imported from google trends: Must be inputted as column vectors
time=Month;
views=Views;

%Initialize parameters with possible values
B=rand(3,1)*100;

%Ode Solver
f=@(t,y) [-B(1)*y(1)*y(2) ; ((B(1)-B(2))*y(1)*y(2))+((B(2)*N-B(3)-(B(2)*y(2))*y(2)))];
trange=[0,length(Month)];
initials=[S0 I0];
[t1, infect]=ode45(f,trange,initials);

[Bnew, Rsdnrm, Rsd, ExFlg, OptmInfo, Lmda, Jmat]=lsqcurvefit(infect(:,2),B,time,views);







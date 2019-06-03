%Parameter Values
alpha=1.62e-4;
beta=1.52e-4;
gamma=3.00e-2;
N=200;
I0=2;
S0=N-I0;

%Ode Solver
f=@(t,y) [-alpha*y(1)*y(2) ; ((alpha-beta)*y(1)*y(2))+((beta*N-gamma-(beta*y(2))*y(2)))];
trange=[0,5000];
initials=[S0 I0];
[time, infect]=ode45(f,trange,initials);

%Plot of Views Per Day
plot(time, infect(:,2),'k-')
xlabel('Time(Day)')
ylabel('Search Volume Index')
title('Search Text = "Blog"')
ylim([0 12]);
xlim([0 5000]);
legend('I(t)')

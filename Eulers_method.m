% Initial concentration of hydrogen peroxide 
h2o2_0 = 2.0 ; 
%Temperature = 333 K 
T=333; 
% Rate constant at T = 333 K 
k =rateconstant(T);  % s^-1 
% Time interval i.e of the 40s 
tspan = [0, 40];  % seconds 
dt = 0.01;  % s 
% Number of time steps 
n = ((tspan(2) - tspan(1)) / dt); 
% Initialize arrays to store concentrations and time 
t = zeros(n + 1, 1); 
h2o2 = zeros(n + 1, 1); 
h2o = zeros(n + 1, 1); 
o2 = zeros(n + 1, 1); 
% Set initial concentration for the reaction 
t(1) = tspan(1); 
h2o2(1) = h2o2_0; 
% Performing Euler's explicit method 
for  i = 1:n 
t(i+1) = t(i) + dt; 
% Calculate the derivatives 
dh2o2_dt = -2 * k * h2o2(i); 
dh2o_dt = 2 * k * h2o2(i); 
do2_dt = k * h2o2(i); 
% Update the concentrations 
h2o2(i+1) = h2o2(i) + dh2o2_dt*dt; 
h2o(i+1) = h2o(i) + dh2o_dt*dt; 
o2(i+1) = o2(i) + do2_dt*dt; 
% Ensuring concentrations of the reactant, as well as the product, don't go below zero 
if  h2o2(i+1) < 0 
h2o2(i+1) = 0; 
end 
if  h2o(i+1) < 0 
h2o(i+1) = 0; 
end 
if  o2(i+1) < 0 
o2(i+1) = 0; 
end 
end 
% Plot the concentrations over time 
figure; 
plot(t, h2o2,  'r-' ,  'LineWidth' , 2); 
hold  on ; 
plot(t, h2o,  'b-' ,  'LineWidth' , 2); 
plot(t, o2,  'g-' ,  'LineWidth' , 2); 
hold  off ; 
xlabel( 'Time (s)' ); 
ylabel( 'Concentration (M)' ); 
legend( 'H2O2' ,  'H2O' ,  'O2' ); 
title( 'Decomposition of Hydrogen Peroxide Simulation' ); 
grid  on ; 
function  k=rateconstant(T) 
Ea=75000; 
A=1.2*(10^11); 
R=8.314; 
k=zeros(size(T,1),1); 
for  i=1:size(T,1) 
k(i)=A*(10^(-Ea/(2.303*R*T(i)))); 
end 
end
clear
clc

% Simulation Parameters
simtime=100; %ms
dt=.01; %ms time step
t=0:dt:simtime; %ms time array

% Problem 1
I(1:length(t))=0; %injection current

% Problem 2
% I(1:length(t))=0;
% I(101:150)=5;

% Problem 3
% I(1:length(t))=5;


% Constants
gK = 36; % mS/cm^2
gNa = 120; % mS/cm^2
gL = 0.3; % mS/cm^2
EK = -12; % mV
ENa = 115; % mV
EL = 10.6; % mV
Cm = 1.0; % uF/cm^2
Vm = 0;

% Equations

% Initial gating variables
alpha_m = 0.1*((25-Vm)/(exp((25-Vm)/10)-1));
beta_m = 4*exp(-Vm/18);
alpha_n = 0.01*((10-Vm)/(exp((10-Vm)/10)-1));
beta_n = 0.125*exp(-Vm/80);
alpha_h = 0.07*exp(-Vm/20);
beta_h = 1/(exp((30-Vm)/10)+1);

% Initial m/n/h values
m(1) = alpha_m/(alpha_m + beta_m);
n(1) = alpha_n/(alpha_n + beta_n);
h(1) = alpha_h/(alpha_h + beta_h);

for i=1:length(t)-1
    
    % Calculating the new alpha and beta values for m/n/h
    alpha_m = 0.1*((25-Vm(i))/(exp((25-Vm(i))/10)-1));
    beta_m = 4*exp(-Vm(i)/18);
    alpha_n = 0.01*((10-Vm(i))/(exp((10-Vm(i))/10)-1));
    beta_n = 0.125*exp(-Vm(i)/80);
    alpha_h = 0.07*exp(-Vm(i)/20);
    beta_h = 1/(exp((30-Vm(i))/10)+1);
    
    % Currents
    INa = (m(i)^3)*gNa*h(i)*(Vm(i) - ENa);
    IK = (n(i)^4)*gK*(Vm(i)-EK);
    IL = gL*(Vm(i)-EL);
    I_ion = I(i) - IK - INa - IL;
    
    % Difference between the current m/n/h value and the next one as well
    % as the membrane voltage and the next voltage value
    dm_dt = alpha_m*(1-m(i)) - beta_m*m(i);
    dn_dt = alpha_n*(1-n(i)) - beta_n*n(i);
    dh_dt = alpha_h*(1-h(i)) - beta_h*h(i);
    dv_dt = I_ion/Cm;

    % Euler method to find the next m/n/h value    
    m(i+1) = m(i) + dt*dm_dt;
    n(i+1) = n(i) + dt*dn_dt;
    h(i+1) = h(i) + dt*dh_dt;
    Vm(i+1) = Vm(i) + dt*dv_dt;
    
end

Vm = Vm-70; % Changing the rest voltage afterward

plot(t,Vm)
xlabel('Time (ms)')
ylabel('Voltage (mV)')
axis([0, 100, -100, 50])
title('Membrane Potential')

% Plotting Conductances
figure
pgK = plot(t,gK*n.^4);
hold on
pgNa = plot(t,gNa*(m.^3).*h,'g');
legend([pgK, pgNa], 'gK', 'gNa')
ylabel('Conductance (mS/cm^2)')
xlabel('Time (ms)')
title('gK and gNa')


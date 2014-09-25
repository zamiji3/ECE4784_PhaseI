clear
clc

% Simulation Parameters
simtime=100; %ms
dt=.04; %ms time step
t=0:dt:simtime; %ms time array


% Constants
gK=36; %mS/cm^2
gNa=120; %mS/cm^2
gL=.3; %mS/cm^2
EK=-12; %mV
ENa=10.6; %mV
EL=10.6; %mV
Vrest=-70; %mV
Cm=1.0; %uF/cm^2
Vm(1)=-70; %mV initial membrane voltage
I=0; %injection current

% Equations

% Gating variables
alpha_m=.1*((25-Vm(1))/(exp((25-Vm(1))/10)-1));
beta_m=4*exp(-Vm(1)/18);
alpha_n=.01*((10-Vm(1))/(exp((10-Vm(1))/10)-1));
beta_n=.125*exp(-Vm(1)/80);
alpha_h=.07*exp(-Vm(1)/20);
beta_h=1./(exp((30-Vm(1))/10)+1);

m(1)=alpha_m./(alpha_m+beta_m);
n(1)=alpha_n./(alpha_n+beta_n);
h(1)=alpha_h./(alpha_h+beta_h);

% Currents
INa=m(1)^3*gNa*h(1)*(Vm(1)-ENa);
IK=n(1)^4*gK*(Vm(1)-EK);
IL=gL*(Vm(1)-EL);
I_ion=I-IK-INa-IL;     

% Derivatives
dVm_dt= I_ion/Cm;
dm_dt=alpha_m*(1-m)-beta_m*m;
dn_dt=alpha_n*(1-n)-beta_n*n;
dh_dt=alpha_h*(1-h)-beta_h*h;

for i=1:length(t)-1
    
    m(i+1)=m(i)+dt.*dm_dt;
    n(i+1)=n(i)+dt.*dn_dt;
    h(i+1)=m(i)+dt.*dh_dt;
    
    alpha_m=.1*((25-Vm(i))/(exp((25-Vm(i))/10)-1));
    beta_m=4*exp(-Vm(i)/18);
    alpha_n=.01*((10-Vm(i))/(exp((10-Vm(i))/10)-1));
    beta_n=.125*exp(-Vm(i)/80);
    alpha_h=.07*exp(-Vm(i)/20);
    beta_h=1./(exp((30-Vm(i))/10)+1);

    % Currents
    INa=m(i)^3*gNa*h(i)*(Vm(i)-ENa);
    IK=n(i)^4*gK*(Vm(i)-EK);
    IL=gL*(Vm(i)-EL);
    I_ion=I-IK-INa-IL;     

    
    Vm(i+1)=Vm(i)+dt*dVm_dt;
end

plot(t,Vm,'g',t,gK,'b')

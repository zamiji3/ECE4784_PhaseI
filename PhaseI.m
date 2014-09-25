% Simulation Parameters
simtime=100; %ms
dt=1; %ms time step
t=0:dt:simtime; %ms time array


% Constants
gK=36; %mS/cm^2
gNa=120; %mS/cm^2
gL=.3; %mS/cm^2
EK=-12; %mV
ENa=10.6; %mV
EL=10.6; %mV
Vrest=-70; %mV

Vm=-70; %mV initial membrane voltage

% Equations

% Gating variables
alpha_m=.1*((25-Vm)/(exp((25-Vm)/10)-1));
beta_m=4*exp(-Vm/18);
alpha_n=.01*((10-Vm)/(exp((10-Vm)/10)-1));
beta_n=.125*exp(-Vm/80);
alpha_h=.07*exp(-Vm/20);
beta_h=1/(exp((30-Vm)/10)+1);

m_o=alpha_m/(alpha_m+beta_m);
n_o=alpha_n/(alpha_n+beta_n);
h_o=alpha_h/(alpha_h+beta_h);

% Currents
INa=m^3*gNa*h*(Vm-ENa);
IK=n^4*gK*(Vm-EK);
IL=gL*(Vm-EL);
I_ion=I-IK-INa-IL;

% Derivatives
dVm_dt= I_ion/Cm;
dm_dt=alpha_m*(1-m)-beta_m*m;
dn_dt=alpha_n*(1-n)-beta_n*n;
dh_dt=alpha_h*(1-h)-beta_h*h;


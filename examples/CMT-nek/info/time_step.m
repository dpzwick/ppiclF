% Computes max time step based on collision scale
rhop = 2500;
dp = 115e-6;
k = 10000;
er = 0.9;
nres = 10;
volp = pi/6*dp^3;
massp = volp*rhop;
rm12 = 1/(1/massp + 1/massp);
dt = sqrt(rm12/k*(log(er)^2+pi^2))/nres

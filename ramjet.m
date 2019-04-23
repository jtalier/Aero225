%% Dan Maguire
%% AEROSP 225
%% Project

clear;
clc;
close all;
format short;

%% Inputs
M1      = 3.00;         % unitless             Mach #
height  = 30E3;         % m                    Cruise Altitude
% Fuel Type: Hydrogen
q_HV    = 120E6;        % J/kg                 Heating Value
MW      = 28.9;         % g/mol                Molecular Weight
MW_fuel = 2;            % g/mol                Molecular Weight
gamma   = 1.4;          % unitless             Specific Heat Ratio
m_dot   = 100;          % kg/s                 Mass Flow Rate
T5      = 1800;         % K                    Exit Combustor Temperature


% Requirements:
% Inlet efficiency > 0.8
% Thrust > 60E3 N
% As small as possible

Ru      = 8.3144598;     % J/(mol*K)           Universal Gas Constant

% Calculation of Gas Constants
R_fuel     = Ru/MW_fuel;        % J/(g*K)
R_fuel     = R_fuel*1000;       % J/(kg*K)
cp_fuel    = gamma/(gamma-1)*R_fuel; % J/(kg*K)
cv_fuel    = 1/(gamma-1)*R_fuel;     % J/(kg*K)

R     = Ru/MW;                  % J/(g*K)
R     = R*1000;                 % J/(kg*K)
cp    = gamma/(gamma-1)*R;      % J/(kg*K)
cv    = 1/(gamma-1)*R;          % J/(kg*K)




%% Initial State
[T1, a1, p1, rho1] = atmoscoesa(height);
[Mrat, Trat, prat, rhorat, arearat] = flowisentropic(gamma, M1);
T01 = T1/Trat;     p01 = p1/prat;   rho01 = rho1/rhorat;
u1 = M1*a1;


%% Inlet
disp('Inlet:');

syms Bsym;

%{
theta = 2;     % deg

M = zeros(1,5);
M(1) = M1;
B = zeros(1,5);
for i = 1:4
    B(i) = fzero(@(B) tand(theta) - 2*cotd(B)* ...
        (M(i)^2*(sind(B))^2 - 1) / (M(i)^2*(gamma + cosd(2*B)) + 2), theta);
    Mn = M(i)*sind(B(i));
    Mn2 = sqrt((Mn^2 + 2/(gamma-1))/((2*gamma/(gamma-1))*Mn^2 - 1));
    M(i+1) = Mn2/sind(B(i)-theta);
end

B
%}

% OS 1
numShocks = 14;
M = zeros(1,numShocks+1);

p = zeros(1,numShocks+1);
p(1) = p1;
p0 = zeros(1,numShocks+1);
p0(1) = p01;
T = zeros(1,numShocks+1);
T(1) = T1;
rho = zeros(1,numShocks+1);
rho(1) = rho1;

M(1) = M1;
B = zeros(1,numShocks);
%theta = zeros(1,numShocks);

theta = 4;     % deg
starting_guess = 60;

for i = 1:numShocks
    if i == 12
        theta = 2
        starting_guess = 60;
    end
    if i == 14
        theta = 1
        starting_guess = 60;
    end
    
    B(i) = fzero(@(B) tand(theta) - 2*cotd(B)* ...
        (M(i)^2*(sind(B))^2 - 1) / (M(i)^2*(gamma + cosd(2*B)) + 2), starting_guess);
    Mn = M(i)*sind(B(i));
    Mn2 = sqrt((Mn^2 + 2/(gamma-1))/((2*gamma/(gamma-1))*Mn^2 - 1));
    M(i+1) = Mn2/sind(B(i)-theta);
    
    prat = 1 + (2*gamma)/(gamma+1)*(Mn^2 - 1);
    p(i+1) = p(i)*prat;
    
    rhorat = (gamma+1)*Mn^2/((gamma-1)*Mn^2 + 2)
    rho(i+1) = rho(i)*rhorat;
    
    Trat = prat/rhorat;
    T(i+1) = T(i)*Trat;
    
    [Mrat, Trat, prat, rhorat, arearat] = flowisentropic(gamma, M(i+1));
    p0(i+1) = p(i+1)/prat;
    
end
 
B
M
p
p0
p0(end)/p0(1)


    

[mach, Trat, prat, rhorat, downstream_mach, p0rat] = flownormalshock(gamma, M(end));
M2 = downstream_mach;
T2 = Trat*T(end);  p2 = prat*p(end);   rho2 = rhorat*rho(end);

[Mrat, Trat, prat, rhorat, arearat] = flowisentropic(gamma, M2);
T02 = T(end)/Trat;      p02 = p(end)/prat;
    
p02/p01


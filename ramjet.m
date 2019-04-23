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
theta = 20;     % deg

fzero(@(B) tand(theta) - 2*cotd(B)* ...
    (M1^2*(sind(B))^2 - 1) / (M1^2*(gamma + cosd(2*B)) + 2), 20)

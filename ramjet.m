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
T5_max  = 1800;         % K                    Exit Combustor Temperature


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

w = 2;  % m, design value, width / depth into page


%% Initial State
[T1, a1, p1, rho1] = atmoscoesa(height);
a1 = sqrt(gamma*R*T1);
[Mrat, Trat, prat, rhorat, arearat] = flowisentropic(gamma, M1);
T01 = T1/Trat;     p01 = p1/prat;   rho01 = rho1/rhorat;
u1 = M1*a1;
A1 = m_dot / rho1 / u1;
h1 = A1 / w;


%% Inlet
disp('Inlet:');


% OS 1
numShocks = 4;
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

theta = 11;     % deg
starting_guess = 30;

for i = 1:numShocks
    if i == 4
        %theta = 11;
        starting_guess = 50;
    end
    
    B(i) = fzero(@(B) tand(theta) - 2*cotd(B)* ...
        (M(i)^2*(sind(B))^2 - 1) / (M(i)^2*(gamma + cosd(2*B)) + 2), ...
        starting_guess);
    Mn = M(i)*sind(B(i));
    Mn2 = sqrt((Mn^2 + 2/(gamma-1))/((2*gamma/(gamma-1))*Mn^2 - 1));
    M(i+1) = Mn2/sind(B(i)-theta);
    
    prat = 1 + (2*gamma)/(gamma+1)*(Mn^2 - 1);
    p(i+1) = p(i)*prat;
    
    rhorat = (gamma+1)*Mn^2/((gamma-1)*Mn^2 + 2);
    rho(i+1) = rho(i)*rhorat;
    
    Trat = prat/rhorat;
    T(i+1) = T(i)*Trat;
    
    [Mrat, Trat, prat, rhorat, arearat] = flowisentropic(gamma, M(i+1));
    p0(i+1) = p(i+1)/prat;
    
end

B
M
%p0(end)/p0(1)


[mach, Trat, prat, rhorat, downstream_mach, p0rat] = ...
    flownormalshock(gamma, M(end));
%% State 3
M3 = downstream_mach;
T3 = Trat*T(end);  p3 = prat*p(end);   rho3 = rhorat*rho(end);

[Mrat, Trat, prat, rhorat, arearat] = flowisentropic(gamma, M3);
T03 = T3/Trat;     p03 = p3/prat;   rho03 = rho3/rhorat;
a3 = sqrt(gamma*R*T3);
u3 = M3*a3;

p03/p01



%% Inlet Geometry
x = 0:0.01:13;
lower_wall_y = x*tand(theta);
shock1_y = x*tand(B(1));

plot(x,lower_wall_y,'k',x,h1*ones(1,length(x)),'k',x,shock1_y);
hold on;
%x*sind(B(1)) = h
hit1_x = h1 / tand(B(1));
% shock 2 m = -tand(B(2)-theta)
% shock 2 point : hit1_x, h
% y - h = -tand(B(2)-theta)*(x-hit1_x)
shock2_y = h1 - tand(B(2)-theta)*(x - hit1_x);
plot(x,shock2_y);

% shock2_y = hit2_y = tand(theta)*hit2_x = h1 -tand(B(2)-theta)*(hit2_x - hit1_x)
hit2_x = (h1+tand(B(2)-theta)*hit1_x) / (tand(theta) + tand(B(2)-theta));

% shock 3 m = sind(B(3))
% shock 3 point : hit2_x, sind(theta)*hit2_x
% y - sind(theta)*hit2_x = sind(B(3))*(x-hit2_x)
shock3_y = tand(theta)*hit2_x +tand(B(3))*(x - hit2_x);

plot(x,shock3_y);


% shock3_y = hit3_y = h1 = sind(theta)*hit2_x +sind(B(3))*(x - hit2_x)
hit3_x = (h1-tand(theta)*hit2_x + tand(B(3))*hit2_x) / ...
    (tand(B(3)));

% shock 4 m = -sind(B(4)-theta)
% shock 4 point : hit3_x, h1
% y - h1 = -sind(B(4)-theta)*(x-hit3_x)
shock4_y = h1 -tand(B(4)-theta)*(x-hit3_x);

plot(x,shock4_y);

% shock4_y = hit4_y = sind(theta)*hit4_x = h1 -sind(B(4))*(hit4_x-hit3_x);
hit4_x = (h1+tand(B(4)-theta)*hit3_x) / (tand(theta) + tand(B(4)-theta));

figure(2);
hold on;
% lower wall
x = 0:0.01:hit4_x;
plot(x,x*tand(theta),'k');  

% upper wall
x = hit1_x:0.01:13;
plot(x,h1*ones(1,length(x)),'k');           

% first shock
x = 0:0.01:hit1_x;
shock1_y = x*tand(B(1));
plot(x,shock1_y); 

% second shock
x = hit1_x:0.01:hit2_x;
shock2_y = h1 -tand(B(2)-theta)*(x - hit1_x);
plot(x,shock2_y);

% third shock
x = hit2_x:0.01:hit3_x;
shock3_y = tand(theta)*hit2_x +tand(B(3))*(x - hit2_x);
plot(x,shock3_y);

% fourth shock
x = hit3_x:0.01:hit4_x;
shock4_y = h1 -tand(B(4)-theta)*(x-hit3_x);
plot(x,shock4_y);

% end of lower wall
x = hit4_x:0.01:13;
plot(x,(tand(theta)*hit4_x)*ones(1,length(x)),'k');


axis([0 13 0 4]);
axis equal;


h3 = h1 - (tand(theta)*hit4_x);
A3 = h3*w;

%% DIFFUSER
A3 = ; %Starting aera of diffuser
h3 = A3/w; %Starting Diffuser Height
A4 = A1; %End area of diffuser
h4 = A4/w; %Diffuser Height

%Find A*
[Mrat, Trat, prat, rhorat, arearat] = flowisentropic(gamma, M3);
a_star = arearat/A3; 

L_diff = 1; %[meters] %Length of diffuser
L = linspace(0,L);
h = (h4-h3/L_diff)*L + h3;




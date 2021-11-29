# This code was written in Octave, https://www.gnu.org/software/octave/index

# This code simulates the motion of solar wind around the Earth.
# Solar wind is modeled as a bunch of noninteracting charged point particles.
# The Earth is modeled as a magnetic dipole and a point mass
# The atmosphere is included through drag force.

close all
clear
clc

%% Physical constants
M = 1; %% Mass of the Earth
Re = 2; ## Earth radius, used to model the atmospheric drag onto particles.
m = 1; ## Mass of the charged particles
q = 0.01; ## Charge of the charged particles
G = 1; ## Newton's gravitational constant
v0 = 0.2; ## magnitude of the initial velocity vector
B0 = 1.2; ## The scale of the Earth's magnetic field

%% Integrator setup
numIter = 1e4; ## Number of iterations
T = 100; ## Length of time simulated.
t = linspace(0, T, numIter); ## time array
h = t(2) - t(1); ## Time step
%% Physics
## Magnetic field is modeled by a simple dipole.
## Earth's gravitational field is modeles as a point source
## Earth's atmospheric drag is modeled as F_drag = -g*exp(-r/Re)*v. The atmosphere is assumed to get exponentionaly thinner and the drag grows linearly with velocity and proportional with atmospheric pressure.

## Atmosphere model
g = 0.22;
## Dipole magnetic field
magnetic_dipole_moment_direction = [0, 0, -1];
B = @(R) B0*(3*sum(magnetic_dipole_moment_direction.*R) * R/norm(R)^2 - magnetic_dipole_moment_direction)/norm(R)^3;

# Equations of motions that include gravity. I think gravity is negligible here.
#F = @(S) [S(2, :); -M*G*S(1, :)/norm(S(1, :))^3 + q /m * cross(S(2, :), B(S(1, :))) - g*exp(-S(1))*S(2, :)]; % a = -M G r/r^3 + q/m (V x B)
# Equations of motion that include only drag and magnetic field
F = @(S) [S(2, :); q /m * cross(S(2, :), B(S(1, :))) - g*exp(-S(1))*S(2, :)]; % a = M G r/r^3 + q/m (V x B)


## The simulation is by making num_initial_conditions particles come as Solar wind.
## The particles are released from a rectangular window
num_initial_conditions = 50;

y0min = -0.24+0.1;
y0max = -0.26+0.1;
z0min = -0.35+0.1;
z0max = -0.34+0.1;

y0min = -1/3;
y0max = +1/3;
z0min = -0.11;
z0max = -0.09;

initial_conditions_data = [y0min + (y0max - y0min)*rand(num_initial_conditions);
                           z0min + (z0max - z0min)*rand(num_initial_conditions);];
plotData = zeros(3, numIter, num_initial_conditions);

for ii = 1 : num_initial_conditions
  S = [1, initial_conditions_data(1, ii), initial_conditions_data(2, ii); -v0, 0, 0];
  plotData(:, 1, ii) = S(1, :);
  for s = 2 : numIter
    k1 = h * F(S);
    k2 = h * F(S + k1);
    k3 = h * F(S + (3*k1 + k2)/8);
    k4 = h * F(S + (8*k1 + 2*k2 + 8*k3)/27);
    k5 = h * F(S + (3*(3*sqrt(21) - 7)*k1 - 8*(7 - sqrt(21))*k2 + 48*(7 - sqrt(21))*k3 - 3*(21-sqrt(21))*k4)/392);
    k6 = h * F(S + (-5*(231 + 51*sqrt(21))*k1 - 40*(7 + sqrt(21))*k2 - 320*sqrt(21)*k3 + 3*(21 + 121*sqrt(21))*k4 + 392*(6 + sqrt(21))*k5)/1960);
    k7 = h * F(S + (15*(22 + 7*sqrt(21))*k1 + 120*k2 + 40*(7*sqrt(21) - 5)*k3 - 63*(3*sqrt(21) - 2)*k4 - 14*(49 + 9*sqrt(21))*k5 + 70*(7 - sqrt(21))*k6)/180);

    t = t + h;
    S = S + (9*k1 + 64*k3 + 49*k5 + 49*k6 + 9*k7)/180;
    plotData(:, s, ii) = S(1, :);
  end
end

if numIter > 4e2
  sampleData = 1 : ceil(numIter/4e2) : numIter;
else
  sampleData = 1 : numIter;
end

figure(1)
plot3(0,0,0, 'r*')
hold on
for ii = 1 : num_initial_conditions
  plot3(plotData(1, sampleData, ii), plotData(2, sampleData, ii), plotData(3, sampleData, ii), 'b-')
  hold on
endfor
xlim([-1/2,1])
ylim([-1/2,1/2])
zlim([-1/2,1/2])
axis square
legend('Dipole magnetic field source','Solar wind particle trajectories')

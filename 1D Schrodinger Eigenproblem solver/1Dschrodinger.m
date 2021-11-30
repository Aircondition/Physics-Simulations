# This code finds the eigenvalues and eigenvectors of 1D TISE problem.

close all
clear
clc

%Initialize the simulation

# define the x axis going from a to b
a = -3;
b = 3;
rezX = 1e4; # spatial resolution of the problem. A higher value generates higher accuracy solutions.
numIter = 20; # number of iterations used in my algorithm FECTP. A higher value generates higher accuracy solutions.
# FECTP means FindEigenvalueClosestToP
x = linspace(a, b, rezX);
dx = x(2) - x(1);

# mass and the Plank's constant
m = 1;
hBar = 1;

U = x.^4 - 5*x.^2-2*x; # Define the potential U(x). This potential has two global minima which should generate some interesting eigenfunctions.
# U = x.^2 SHO
# The wavefunction is seen as a vector with rezX components. The hamiltonian acts via a second derivative operator which is approximated as a second difference matrix in this vector space.

SDM = sparse(rezX); # Second Difference Matrix is -2 on the diagonal and 1 on the first offdiagonal
p = 1 : rezX; # I like using array indexing in MATLAB because for loops should be avoided when possible
SDM(p + (p-1)*rezX) = -2;
SDM((p(1:end-1)+1) + (p(1:end-1)-1)*rezX) = 1;
SDM((p(2:end)-1) + (p(2:end)-1)*rezX) = 1;

temp = sparse(rezX);
temp(p + (p-1)*rezX) = U(p); # Diagonal elements of the hamiltonian contain the potential in this representation.
H = reshape(-0.5 * hBar^2/2/m * SDM/dx^2 + temp, [rezX, rezX]);

%% Search for Eigenvalues and Eigenvectors. Define the range of energies to be searched. From Emin to Emax with resolution numE.
Emax = 10; 
Emin = min(min(U)); % The eigenvalues must be greater than the minimum of the potential so that is a good place to start the search
numE = 5e2;
# Defining the search space of energy levels L
L = linspace(Emin, Emax, numE);
dE = L(2) - L(1);
# Defining an array that will contain the eigenvalues closest to the one in the search space.
data = zeros(1, length(L));
for s = 1 : length(L)
  # pick a candidate eigenvalue, L(s).
  # execute FECTP on this hamiltonian for numIter number of times.
  # The result is an array of length rezX+1 whose last element is the eigenvalue and first rezX elements are the eigenfunction.
  # This eigenvalue is one closest to L(s)
    FECTP_output = FECTP(H, L(s), numIter); 
    data(s) = FECTP_output(end); # put the closest eigenvalue into data array
end

# if you plot data array you will see it contain plateaus. The height of the plateau is (very likely) the correct eigenvalue 
# The rest of the code detects these plateaus and assigns them the mean value of all members of the plateau.
# I define plateau as a subarray of data that contains succesive elements of almost equal value. By almost equal I mean they are within dE from each other.

top = abs(diff(data)) < dE;
top(1) = 0;
blj = find(top == 0); # find all of the places the plateau begins and ends

# This part of the code finds the plateaus' energies more precisely. Each energy is labeled with a quantum number n
n = 0;
marker1 = 1;
marker2 = 1;

for s = 2 : length(blj)
        marker1 = marker2;
        marker2 = blj(s);
        if marker2 - marker1 > 5 # if there are more that 5 succesive plateau points, declare that to be an actual valid energy level
            n = n + 1;
            marker1 = marker1+1;
            EnergyLevel(n) = sum(data(marker1:marker2))/length(marker1:marker2); # average eigenvalue in the plateau
            AbsoluteError(n) = std(abs(EnergyLevel(n) - data(marker1:marker2))); # here I estimate the error in found energy by taking the standard deviation of the plateau energies
        end
end
figure(1)
subplot(1,2,1)
plot(L, data, 'r-')
ylim([-1+min(U), 1+max(EnergyLevel)])
legend('Output of FECTP algorithm')
xlabel('guessed E')
ylabel('output of FECTP(E)')
subplot(1,2,2)
plot(L, data, 'r-')
ylim([-1+min(U), 1+max(EnergyLevel)])
hold on
for ii = 1 : length(EnergyLevel)
    plot(L, EnergyLevel(ii)+0*L, 'k-', L, EnergyLevel(ii) + AbsoluteError(ii)+0*L, 'b-', L, EnergyLevel(ii) - AbsoluteError(ii)+0*L, 'b-')
    legend('Output of FECTP algorithm', 'Average', 'Error bars')
    hold on
end
xlabel('guessed E')
ylabel('output of FECTP(E)')
title('Plateaus of energy')
%% This part of the code refines the solutions found.
figure(2)
numIter = 1e2;
plotScale = 0.5*max(abs(diff(EnergyLevel)));
for s = 1 : length(EnergyLevel)
    kek = FECTP(H, EnergyLevel(s), numIter);
    plot(x, kek(end) + plotScale*kek(1:end-1)/(max(kek(1:end-1)) - min(kek(1:end-1))), "linewidth", 1.2, x, kek(end) + 0*x)
    hold on
end
plot(x, U, "linewidth", 1.2, 'k-')
ylim([-1+min(U), 1+max(EnergyLevel)])
xlim([a, b])
title('Eigenvalues and Eigenvectors of the Hamiltonian H = p^2/2m + V(x)')

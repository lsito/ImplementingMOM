%% Implementing MOM for a lambda/2 dipole antenna

%% Cleaning 
clear;
clc;
close all; 

%% Eelctroamgnetic parameters
c_0 = 3e8;
freq = 300e6;
%freq=10e6;
omega = 2*pi*freq;
lambda = c_0/freq;
kappa = 2*pi/lambda;
eps_0 = 8.854e-12;
mi_0 = 4*pi*1e-7;

%% Alimetation
Vg = 1;

%% Antenna's parameters
a = lambda/100;
%a=lambda/20;
length = lambda/2;
l = length/2;

%% Defining delta_z and the matching points
delta_z = lambda/30;
z_M = linspace(-l, l - delta_z, 2*l/delta_z); %Vector containing only the left element of each interval
z_MP = z_M+(delta_z/2); %Vector containing the matching points

%% Weights and nodes
w_1 = 0.5;
w_2 = 0.5;
x_1 = (-1/sqrt(3)+1)*0.5;
x_2 = (1/sqrt(3)+1)*0.5;

%% Filling the martix A 
n = size(z_M); %Number of columns of the matrix, it gives a 2 elements vector
A = zeros(n(2)); %Creating a matrix of zeros !n(2) is the second element of n

for m = 1:n(2)  
    for k = 1:n(2)
        if k == m
           A(m,m)= (1/(4*pi))*(mi_0/(4*pi))*log((sqrt(1+((4*a^2)/delta_z^2))+1)/(sqrt(1+(4*a^2)/delta_z^2)-1))-((1i*kappa*delta_z)/(4*pi));
        else
            A(k,m) = (delta_z)*(mi_0/(4*pi))*((w_1*integranda(a, z_M(m), z_MP(k), delta_z, x_1, kappa))+(w_2*integranda(a, z_M(m), z_MP(k), delta_z, x_2, kappa)));
        end
    end
end

% for m=1:n(2)
%     for k=1:n(2)
%         if k==m
%       
%           A(m,m)= (1/(4*pi))*log((sqrt(1+((4*a^2)/delta_z^2))+1)/(sqrt(1+(4*a^2)/delta_z^2)-1))-((1i*kappa*delta_z)/(4*pi));
%         else
%            f=@(z3) (delta_z)* mi/(4*pi)*exp(-i*k*sqrt(a.^2+(z_mp(m)-(z3*delta_z+zm(n)))^2))/(sqrt(a^2+(z_mp(m)-(z3*delta_z+zm(n)))^2));
%           
%           A(m,n)=w1*f(x1)+w2*f(x2);
%     
%         end
%     end
    


%% Known terms vector

%Implementing the two addenda
c = cos(kappa*z_MP);
s = -1i*omega*eps_0*mi_0*Vg/(2*kappa)*sin(kappa*abs(z_MP));

%Implementing c (constant)
A_inv = inv(A);
% A_inv = A/eye(n(2)); Other way to invert the matrix

%Base column vector u
u = zeros(n(2),1); 
u(1) = 1;

%Constant C
C = -(u'*A_inv*s')/(u'*A_inv*c');

%Currents' vector
I = C*A_inv*c'+A_inv*s';

%% Plotting
I_real = abs(I);
xlabel('z');
plot(I_real) 


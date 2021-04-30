clear;
clc;
close all;

%% Descrizione
% Discretizzo l'input (I(z')) dell'integrale con delle finestre rettangolari
% Discretizzo l'output (z) dell'integrale col point matching (scelgo il centro degli intervalli)

% Il MOM trasforma un'equazione integrale in un sistema di equazioni lineari.
% Le incognite sono i pesi (In) delle funzioni di base con cui esprimo l'incognita (I): AI = T.

% Risolvo l'integrale col metodo di quadratura gaussiana.
% Se k = m, nella matrice A, cioè negli integrandi, l'integrale vale:
% mi_0/(4*pi)*(1/(4*pi)*log((sqrt(1 + 4*a^2/(delta_z)^2) + 1)/(sqrt(1 +
% 4*a^2/(delta_z)^2) - 1)) - 1i*(kappa*delta_z)/(4*pi)).

%% Dati
epsilon_0 = 8.854e-12;
mi_0 = 4*pi*1e-7;
freq = 300e6;
omega = 2*pi*freq;
c = 3e8;
lambda = c/freq;
l = lambda/4;
a = lambda/100;
k = 2*pi/lambda;
Vg = 1;

%% Calcolo le z e le z_c
delta_z = lambda/30;
N = 2*l/delta_z;

z = linspace(-l, l-delta_z, N);
z_m = z + delta_z/2;

%% Fill di A
A = zeros(N);
w1 = 0.5;
w2 = w1;
x1 = (-1/sqrt(3)+1)*0.5;
x2 = (1/sqrt(3)+1)*0.5;
for m = 1:N
    for n = 1:N
        if n == m
            A(m, n) = (mi_0/(4*pi))*log((sqrt(1+((4*a^2)/delta_z^2))+1)/(sqrt(1+(4*a^2)/delta_z^2)-1))-((1i*k*delta_z)/(4*pi));
        else
            A(m, n) = mi_0/(4*pi)*((w1*delta_z*integranda(delta_z, k, a, z_m(n), z(m), x1)) + w2*delta_z*(integranda(delta_z,k,a,z_m(n),z(m),x2)));
        end
    end
end

%% Calcolo il termine noto
u = zeros(1, N);    
u(1) = 1;
u(N) = 1;

c = (cos(k*z_m))';
s = (-1i*((2*pi)*freq)*epsilon_0*Vg)/(2*k)*sin(k*abs(z))';
C = -(u*inv(A)*s)/(u*inv(A)*c);
I = inv(A)*C*c + inv(A)*s;
plot(z/lambda, abs(I));
xlabel('$l(\lambda)$', 'Interpreter', 'latex');
ylabel('$|I(l)|$', 'Interpreter', 'latex');

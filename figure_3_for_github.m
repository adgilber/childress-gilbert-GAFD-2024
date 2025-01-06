% script to  solve 3-equation level 1 model;    s = sigma

% equations for [g, h, th](s)

clear all
close all

level = 1  ;

beta  = 4      ;
C     = 0.5 ;

% following loop can be removed (with the `end' at the bottom of the script

ictype =  0 ;   % 0 is given, 1 is a branch, 2 is b branch, for small sigma expansion

gamma = 1/beta ;
mu    = (beta-1)*gamma + 1 ;

dsigma = 0.0001 ;
sigmax = 20 ;
sigma0 = 0.01   ;

%% standard initial conditions

g0  =  1 ;
h0  =  0 ;
th0 =  0 ;

disp('original initial values')
disp([g0, h0, th0])

D =   1-C ;
E =  (C-1) * ( (2*beta-1)*gamma + 1) - gamma ;
F = - beta*gamma*mu*C ;

h1a = 0.5 ./ D .* ( - E + sqrt( E.^2 - 4 * F .* D ) ) ;
h1b = 0.5 ./ D .* ( - E - sqrt( E.^2 - 4 * F .* D ) ) ;

th1a  = beta*gamma - h1a ;
th1b  = beta*gamma - h1b ;

g2a = 0.5 * (h1a-mu) * th1a ;
g2b = 0.5 * (h1b-mu) * th1b ;

disp('original initial values')
disp([g0, h0, th0])

%% or early-time expansions


if ictype == 1

    g0  =  1   + g2a  * sigma0^2  ;
    h0  =        h1a  * sigma0    ;
    th0 =        th1a * sigma0    ;

    disp('updates, branch a')
    disp([g2a, h1a, th1a])

    disp('updated initial values, branch a')
    disp([g0, h0, th0])

elseif ictype == 2

    g0  =  1   + g2b  * sigma0^2  ;
    h0  =        h1b  * sigma0    ;
    th0 =        th1b * sigma0    ;

    disp('updates, branch b')
    disp([g2b, h1b, th1b])

    disp('updated initial values, branch b')
    disp([g0, h0, th0])

end

% make real if necessary

g0  =  real(g0)  ;
h0  =  real(h0)  ;
th0 =  real(th0) ;

%% start time-stepping

nsteps = round(sigmax / dsigma  - 1)  ;

g    = g0  ;
h    = h0  ;
th   = th0 ;

s    = sigma0 ;
zeta = 0.0    ;
x    = 0.0    ;
y    = 0.0    ;

gstore (1) = g  ;
hstore (1) = h  ;
thstore(1) = th ;

sstore (1) = s  ;
zstore (1) = zeta ;
xstore (1) = x ;
ystore (1) = y ;

condAstore(1) = 1.0 ;
normvstore(1) = 1.0 ;

for n=1:nsteps+1

    %% compile matrix

    A21 = - C * g^(1-beta)    ;

    A31 = g^(-beta) ;

    A = [ mu*s, -1/beta*g^(1-beta),  -1/beta*g^(2-beta)       ; ...
        A21,         mu*s,             g^(1-beta)*h           ; ...
        A31,           0,              mu*s-g^(-beta)*h     ] ;

    b = [ -gamma*g ; -gamma*h ; 0 ] ;

    %% solve matrix problem and step along

    v = A\b ;

    g  = g  + dsigma*v(1) ;
    h  = h  + dsigma*v(2) ;
    th = th + dsigma*v(3) ;

    s  = s  + dsigma      ;

    k  = g^(-beta) * v(3) ;

    dzeta = g^beta  * dsigma ;
    zeta  = zeta + dzeta  ;
    x     = x    + dzeta*cos(th) ;
    y     = y    + dzeta*sin(th) ;

    kstore(n) = k ;

    if n < nsteps+1

        gstore(n+1)  = g  ;
        hstore(n+1)  = h  ;
        thstore(n+1) = th ;

        sstore(n+1)  = s ;
        zstore(n+1)  = zeta ;
        xstore(n+1)  = x ;
        ystore(n+1)  = y ;

        condAstore(n+1) = cond(A) ;
        normvstore(n+1) = norm(v) ; 

    end
end

%% plotting

close all

figure(1)
plot(sstore, gstore, 'b')
hold on 
plot(sstore, -5*hstore, 'r')
xlabel('sigma')
ylabel('g, -5h')
title ('g(sigma), -5h(sigma)')
xlim([0,20])
 
figure(2)
plot(sstore, gstore, 'b')
hold on 
plot(sstore, -5*hstore, 'r')
xlabel('sigma')
ylabel('g, -5h')
title ('g(sigma), -5h(sigma)')
xlim([0,1.0])
 

figure(3)
plot(xstore, ystore, 'b', -xstore, ystore, 'r')
hold on

xlabel('x')
ylabel('y')
title ('hairpin shape')
ylim([0.0,0.2])

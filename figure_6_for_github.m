% script to  solve level N model;    s = sigma

% v25 equations (8.20-8.23)  [g, h, th, Q2, G2, Q3, G3, Q4, G4, ...](s)

clear all
%close all

Nmax   = 4 ;
test4  = 0 ;    %set to 1 to test for Nmax = 4 against previous hard-coded version
level  = Nmax  ;
indmax = 2*Nmax+1 ;

g0  =  1   ;
h0  =  0   ;
h1  =  0   ;
th0 =  0   ;

beta  = 4    ;
gamma = 1/beta ;
C     = 3.35 ;

dsigma = 0.0001 ;
sigmax = 20.0   ;
sigma0 = 0.001  ;

nsteps = round(sigmax / dsigma)  - 1  ;
mu     = (beta-1) * gamma + 1  ;   % this is mu

g    = g0  ;
h    = h0  ;
th   = th0 ;
s    = sigma0 ;
zeta = 0.0 ;
x    = 0.0 ;
y    = 0.0 ;

% pre-allocate all the storage we need...

gstore   = zeros(1,nsteps+1) ;
hstore   = zeros(1,nsteps+1) ;
thstore  = zeros(1,nsteps+1) ;
kstore   = zeros(1,nsteps+1) ; 
sstore   = zeros(1,nsteps+1) ;
zstore   = zeros(1,nsteps+1) ;
xstore   = zeros(1,nsteps+1) ;
ystore   = zeros(1,nsteps+1) ;

Q      = zeros(Nmax,1)    ;
G      = zeros(Nmax,1)    ;
Qstore = zeros(Nmax,nsteps+1)    ;
Gstore = zeros(Nmax,nsteps+1)    ;

% store initial conditions

gstore (1) = g  ;
hstore (1) = h  ;
thstore(1) = th ;
sstore (1) = s  ;
zstore (1) = zeta ;
xstore (1) = x ;
ystore (1) = y ;

for j = 1:Nmax

    Q(j)        = 0.0 ;
    G(j)        = 0.0 ;
    Qstore(j,1) = 0.0 ;
    Gstore(j,1) = 0.0 ;
    
end

%% set up any special initial conditions

% G(2)        =  0.8 ; 
% G(3)        =  0.0 ; 
% G(4)        =  0.223 ;
% Gstore(2,1) = G(2) ; 
% Gstore(3,1) = G(3) ; 
% Gstore(4,1) = G(4) ; 


%% start the actual time-stepping

for n=1:nsteps+1

    A = zeros(indmax,indmax) ;
    b = zeros(indmax,1) ;

    % calculate A21:

    A21 = - C*g^(1-beta) - 2*h0*beta*mu*s*g^(beta-1) + h0^2*beta*g^(beta-1) ;

    if Nmax > 1

        for N = 1:Nmax-1

            A21 = A21 - 0.5*C*((2*N+2)-(2*N+1)*beta)*g^((2*N+1)-(2*N+2)*beta)*G(N+1) ;

        end

    end

    % fill in the 3x3 piece

    A(1,1) = mu*s ;
    A(1,2) = -1.0/beta*g^(1-beta) ;
    A(1,3) = -1.0/beta*g^(2-beta) ;

    A(2,1) = A21  ;
    A(2,2) = mu*s ;
    A(2,3) = g^(1-beta)*h ;

    A(3,1) = g^(-beta) ;
    A(3,3) = mu*s-g^(-beta)*h ;

    b(1) = -gamma*g ;
    b(2) = -gamma*h+h0*(beta+1)*gamma*g^beta ;
    b(3) = -h0*beta*gamma ;

    % fill in the rest of row 2:

    if Nmax > 1

        for N = 1:Nmax-1

            coeff = -0.5*C * g^((2*N+2) - beta*(2*N+2)) ;  % of d_sigma G_{N+1}
            A(2,2*N+3) = coeff ;

        end

    end

    % fill in the Q_2 line

    if Nmax > 1

        A(4,1) = 0.25*C*beta* g^(beta-1) ;
        A(4,4) = mu*s - h0 ;

    end

    % fill in the G_2 line

    if Nmax > 1

        A(5,1) = ( -(mu*s-h0)*G(2) + Q(2) ) * 2*beta*g^(-1) ;
        A(5,4) = -1.0 ;
        A(5,5) = mu*s - h0 ;

    end

    % fill in the Q_j lines

    if Nmax > 2

        for N = 3:Nmax

            row = 2*N ;

            A(row,1)     = ((2-N)*(mu*s-h0) * Q(N) + 0.5*C/N * G(N-1))*beta*g^(-1) ;
            A(row,row)   = mu*s-h0 ;

            for j = 2:N-1

                coeff = (N+1-j)*(2*j-N)/N * g^(-beta) * Q(N+1-j) ; % of d_sigma Q_j
                A(row,2*j) = coeff ;

            end

        end

    end

    % fill in the G_j lines

    if Nmax > 2

        for N = 3:Nmax

            row = 2*N+1;

            A(row,1)     = (-(mu*s-h0)*G(N) + Q(N))*N*beta*g^(-1) ;
            A(row,row-1) = -1.0 ;
            A(row,row)   = mu*s-h0 ;

            for j = 2:N-1

                coeff = - (N+1-j) * g^(-beta) * G(N+1-j) ; % of d_sigma Q_j
                A(row,2*j) = coeff ;

                coeff = + (N+1-j) * g^(-beta) * Q(N+1-j) ; % of d_sigma G_j
                A(row,2*j+1) = coeff ;

            end

        end

    end

    format short
    % if mod(n,1000) == 0 
    % disp(s)
    % disp(cond(A))
    % disp(rcond(A))
    % end

    if rcond(A) < 1e-12

    disp('stop as A is becoming singular')
    break

end


    % fill in b

    if Nmax > 1

        for N = 2:Nmax
            b(2*N)   = -(N*beta-2*N+3) * gamma*Q(N) ;
            b(2*N+1) = -(N-1)*(beta-2) * gamma*G(N) ;
        end

    end

    %% test against level 4 calculations

    if Nmax == 4 & test4 == 1

        G2 = G(2) ; G3 = G(3) ; G4 = G(4) ; Q2 = Q(2) ; Q3 = Q(3) ; Q4 = Q(4) ;

        f = - C*g^(1-beta) - 2*h0*beta*mu*s*g^(beta-1) + h0^2*beta*g^(beta-1) ;

        A21 = f - 0.5*C*(4-3*beta)*g^(-4*beta+3)*G2 ...  % level 2 update
                - 0.5*C*(6-5*beta)*g^(-6*beta+5)*G3 ...  % level 3 update
                - 0.5*C*(8-7*beta)*g^(-8*beta+7)*G4 ;    % level 4 update

        A51 = 2*(-(mu*s-h0)*G2 + Q2)*beta*g^(-1) ;

        A71 = 3*(-(mu*s-h0)*G3 + Q3)*beta*g^(-1) ;

        A91 = 4*(-(mu*s-h0)*G4 + Q4)*beta*g^(-1) ;

        A41 = (1.0/4.0)*C*beta*g^(beta-1) ;

        A61 = (1.0/6.0)*C*beta*G2*g^(-1) -   (mu*s-h0)*beta*Q3*g^(-1) ;

        A81 = (1.0/8.0)*C*beta*G3*g^(-1) - 2*(mu*s-h0)*beta*Q4*g^(-1) ;

        Atest = [ mu*s, -1.0/beta*g^(1-beta), -1.0/beta*g^(2-beta), 0 ,      0,            0,         0,            0,        0             ; ...
            A21,            mu*s,             g^(1-beta)*h,  0 ,     -0.5*C*g^(-4*beta+4), 0, -0.5*C*g^(-6*beta+6), 0, -0.5*C*g^(-8*beta+8) ; ...
            g^(-beta),           0,       mu*s-g^(-beta)*h,  0 ,             0,            0,         0,            0,        0             ; ...
            A41,             0,                 0,        mu*s-h0,           0,            0,         0,            0,        0             ; ...
            A51,             0,                 0,         -1.0,          mu*s-h0,         0,         0,            0,        0             ; ...
            A61,             0,                 0, (2.0/3.0)*g^(-beta)*Q2,   0,         mu*s-h0,      0,            0,        0             ; ...
            A71,             0,                 0,  -2*g^(-beta)*G2, 2*g^(-beta)*Q2,     -1.0,     mu*s-h0          0,        0             ; ...
            A81,             0,                 0,           0 ,             0,         g^(-beta)*Q2, 0,         mu*s-h0,     0             ; ...
            A91,             0,                 0,  -3*g^(-beta)*G3, 3*g^(-beta)*Q3, -2*g^(-beta)*G2, 2*g^(-beta)*Q2, -1,   mu*s-h0         ] ;

        btest = [ -gamma*g ; -gamma*h+h0*(beta+1)*gamma*g^beta ; -h0*beta*gamma ; ...
            -(2*beta-1)*gamma*Q2 ;   -(beta-2)*gamma*G2 ; ...
            -(3*beta-3)*gamma*Q3 ; -2*(beta-2)*gamma*G3 ; ...
            -(4*beta-5)*gamma*Q4 ; -3*(beta-2)*gamma*G4 ] ;

        Aerror = max(max(abs(A-Atest))) ; 
        berror = max(abs(b-btest)) ;

        disp([s, Aerror, berror])

    end

    %% actually solve the matrix problem

    v = A\b ;

    g  = g  + dsigma*v(1) ;
    h  = h  + dsigma*v(2) ;
    th = th + dsigma*v(3) ;

    for N = 2:Nmax
        Q(N) = Q(N) + dsigma*v(2*N)  ;
        G(N) = G(N) + dsigma*v(2*N+1);
    end

    s  = s  + dsigma      ;

    k  = g^(-beta) * v(3) ;

    dzeta = g^beta  * dsigma ;
    zeta  = zeta + dzeta  ;
    x     = x    + dzeta*cos(th) ;
    y     = y    + dzeta*sin(th) ;

    kstore(n) = k ;

    if n < nsteps+1

        gstore (n+1)  = g  ;
        hstore (n+1)  = h  ;
        thstore(n+1)  = th ;
        sstore (n+1)  = s ;
        zstore (n+1)  = zeta ;
        xstore (n+1)  = x ;
        ystore (n+1)  = y ;

        for N = 2:Nmax
            Qstore(N,n+1) = Q(N) ;
            Gstore(N,n+1) = G(N) ;
        end

    end
end

%% 

close all

if level == 1 
    col = 'b' ;
elseif level == 2
    col = 'r' ;
elseif level == 3
    col = 'm' ;
elseif level >= 4
    col = 'k' ;
end

for N = 2:Nmax



    figure(1)
    plot(sstore, Gstore(N,:), col)
    hold on
    xlabel('sigma')
    ylabel('G_N')
    title ('G_N(sigma)')
    xlim([0,20])

    figure(2)
    plot(sstore, Qstore(N,:), col)
    hold on
    xlabel('sigma')
    ylabel('Q_N')
    title ('Q_N(sigma)')
    xlim([0,20])

    disp('j, maximum values of Q_j and G_j')
    
    disp([N, max(Qstore(N,:)), max(Gstore(N,:)) ])

end


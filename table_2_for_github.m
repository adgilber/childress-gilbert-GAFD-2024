
clear all
close all

% Storing eps, Amax, ommin, E1, E2, E, Gamma, P0

% P0 values read from figure 11 of ChGi18 - rather inaccurate but okay for our purposes...

Sdata = ...
    [ 0,    37.11, 1,     59.68,  29.84, 89.52,  36.96, 1.7/38 ; ...
      0.05, 40.82, 0.916, 65.06,  31.60, 96.66,  39.00, 1.7/41 ; ...
      0.1,  45.36, 0.830, 71.48,  33.59, 105.07, 41.40, 1.7/45 ; ...
      0.15, 50.89, 0.745, 79.22,  35.89, 115.11, 44.22, 1.7/51 ; ...
      0.2,  57.78, 0.659, 88.68,  38.56, 127.23, 47.56, 1.7/58 ; ...
      0.25, 66.52, 0.575, 100.48, 41.71, 142.19, 51.55, 1.7/66 ; ...
      0.3,  77.82, 0.493, 115.41, 45.46, 160.87, 56.38, 1.7/78 ; ...
      0.35, 92.75, 0.413, 134.77, 49.99, 184.77, 62.56, 1.7/93 ; ...
     0.375,102.25, 0.374, 146.74, 52.79, 199.54, 65.66, 1.7/102 ; ...
      0.4, 113.63, 0.336, 162.06, 56.00, 218.06, 70.06, 1.7/114 ] ; 

% store data in meaningful arrays 

eps     = Sdata(:,1)' ;
Amaxst  = Sdata(:,2)' ;
omminst = Sdata(:,3)' ;
E1st    = Sdata(:,4)' ;
E2st    = Sdata(:,5)' ;
Est     = Sdata(:,6)' ;
Gammast = Sdata(:,7)' ;
P0st    = Sdata(:,8)' ;

% store values for eps = 0 

Amaxst0  = Amaxst(1)  ;    % this is also K^2 
Est0     = Est(1)     ; 
P0st0    = P0st(1)    ;    % this is also C/2K^2
Gammast0 = Gammast(1) ; 

% calculate quantities as functions of eps and plot

Ucurly = 1./sqrt(omminst) .* (Est0 ./ Est).^0.25 ;               % see (D21) 
Acurly = omminst .* (Amaxst ./ Amaxst0) .* sqrt(Est0 ./ Est) ;   % see (D20) 
Pcurly = 1./omminst .* (P0st ./ P0st0) ;                         % see (D18)
Z      = 1./Acurly .* (1./omminst - 1.0) ;                       % see (D23) 

% recalculate P0 from other data 

P0stnew     = 1./eps .* (1.0-omminst) ./ Amaxst ; 
psimax      = 1./eps .* (1.0-omminst) ;        % interesting to compare with data in the FDR paper 
P0stnew(1)  = P0st0         ;                  % fix up initial case
P0st0new    = P0stnew(1)    ;     
Pcurlynew   = 1./omminst .* (P0stnew ./ P0st0new) ;                         % see (D18)

% calculate scales

K = sqrt(Amaxst(1)) ; 
T = (1.0/K) * omminst ;
L = sqrt(omminst) / K .* (Est0 ./ Est).^0.25 ; 

disp('values of Z')
disp(Z)
disp('values of curly-U')
disp(Ucurly)
disp('values of curly-A')
disp(Acurly)


 
%*******************************************************************************
%Cette fonction algo_SQP implemente l'algorithme SQP
%                                                                              *
%                                                                              *
%                                                                              *
%*******************************************************************************
%                                                    %**************************
%                                                    % PARAMETRES EN ENTREE    *
%                                                    %**************************
%
%            une_f                % une fonction dont on cherche un minimum    %
%
%            des_c                % les contraintes des probleme               %
%
%            un_x0                % un incontournable point initial  
%
%           un_lambda0            % une valeur intiale pour les
%           multiplicateurs de Lagrange.
%
%            un_gf                % une fonction qui code le gradient de une_f %
%
%
%            jac_des_c            % une fonction qui code la jacobienne de des_c %
%
%            un_hf                % une fonction qui code le Hessien de  une_f
%            h_des_c                % une fonction qui code les Hessiens de des_c 
%
%            un_nit_max           % nombre maximum d'iterations autorisees     %
%                                 % risees                                     %
%            une_tol_x            % seuil de stationnarite des x_k             %
%            une_tol_f            % seuil de stationnarite des f_k             %
%            une_tol_g            % seuil validant EULER    en x_k             %
%                                   le sous-probl�me:
%                                   (1 exact, 2 pas de Cauchy, 3 CG tronqu�)   %
%            varargin             % la fonction du noyau MATLAB qui collecte   %
%                                 % les parametres additionnels necessaires au %
%                                 % calcul de une_f, donc un_gf eventuellement %
%                                 % donc une_hf eventuellement                 %
%
%                                                    %**************************
%                                                    % PARAMETRES EN SORTIE    *
%                                                    %**************************
%
%            x_opt                % la solution proposee par trust_region      %
%            f_opt                % une_f (x_opt, varargin{:})                 %
%            g_opt                % un_gf (x_opt, varargin{:})                 %
%            fin                  % la cause de l'arret de l'algorithme        %
%            nit                  % le nombre iterations de trust_region       %
%            f_count              % le nombre d'evaluations de une_f           %
%            g_count              % le nombre d'evaluations de un_gf
%
% Responsable: Y. Diouane (youssef.diouane@isae.fr) --
% (C) Institut Sup�rieur de l'A�ronautique et de l'Espace (ISAE-SUPAERO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [                                                             ...
    x_opt  ,                                                           ...
    f_opt  ,                                                           ...
    fin    ,                                                           ...
    nit    ,                                                            ...
    x_seq]                                                                  ...
    = algo_SQP(                                     ...
    une_f          ,                                ...
    des_c          ,                                ...
    un_x0          ,                                ...
    un_lambda0     ,                                ...
    un_gf          ,                                ...
    jac_des_c      ,                                ...
    un_hf          ,                                ...
    h_des_c        ,                                ...
    un_nit_max     ,                                ...
    une_tol_x      ,                                ...
    une_tol_f      ,                                ...
    une_tol_g                                       ...
    )

%**********************
%    CODE             *
%**********************

global f_count  ;                  % nombre     d'evaluations de la
% fonction a minimiser, sans pitie

global g_count  ;                  % nombre     d'evaluations du
% gradient g_f        , sans pitie
global h_count  ;                  % nombre     d'evaluations de la
% Hessienne h_f       , sans pitie
global c_count  ;                  % nombre     d'evaluations des
% constraintes a minimiser, sans pitie

global jc_count  ;                  % nombre     d'evaluations du
% gradient jac_c        , sans pitie
global hc_count  ;                  % nombre     d'evaluations de la
% Hessienne h_c       , sans pitie

f_count    = 0                                                             ;
g_count    = 0                                                             ;
h_count    = 0                                                             ;

c_count    = 0                                                             ;
jc_count   = 0                                                             ;
hc_count   = 0                                                             ;

% tempx sera desormais    %
% vecteur colonne         %
tempx           = un_x0'                                                   ;
templambda      = un_lambda0'                                              ;

n               = length(tempx   )                                         ;
p               = length(templambda   )                                    ;

% Evaluer notre Lagrangien
fdex            = feval(une_f,tempx)                                       ;
cdex            = feval(des_c,tempx)                                       ;

% Evaluer le gradient de notre Lagrangian
gfdex           = feval(un_gf,tempx)                                       ;
jcdex           = feval(jac_des_c,tempx)                                   ;


% Evaluer le Hessien de notre Lagrangian
hfdex           = feval(un_hf,tempx)                                       ;
hcdex           = feval(h_des_c,tempx)                                     ;

k               = 0                                                        ;
fin             = 0                                                        ;

Pzeros = zeros(p);
x_seq = tempx';
    
while(fin==0)
    
    fdex = feval(une_f,tempx);
    cdex = feval(des_c,tempx);
    
    gfdex = feval(un_gf,tempx);
    jcdex = feval(jac_des_c,tempx);
    
    hfdex = feval(un_hf,tempx);
    hcdex = feval(h_des_c,tempx);

    H_L = hfdex + hcdex*templambda;
    M = [H_L, jcdex'; jcdex, Pzeros];
    B = [-gfdex; -cdex];
    X = M\B;
    dk = X(1:n);
    
    tempx = tempx + dk;
    templambda = X(n+1:end);
    gldex = gfdex + jcdex'*templambda;

    x_seq = [x_seq; tempx'];
    
    
    if (norm(dk) < une_tol_x) && (norm(feval(des_c,tempx)) < 1e-6) && (norm(gldex) < 0.01*une_tol_g)
        fin = 1;
        % Insert patience criteria
    end     
    
    
    if(k==un_nit_max)
        fin = 3;        
    end
    
    k = k + 1;
end
x_opt    =   tempx                                                         ;
f_opt    =    fdex                                                         ;
nit      =   k                                                             ;

%*******************************************************************************
%Cette fonction algo_SQP_BFGS implemente l'algorithme SQP (sans derivee)
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
%            H0                   % une initialisation pour la matrice H_k     %
%
%            une_tol_h            % une tol. pertubation pour approche les
%            derivees (via diffenreces finies)
%
%            un_nit_max           % nombre maximum d'iterations autorisees     %
%                                 % risees                                     %
%            une_tol_x            % seuil de stationnarite des x_k             %
%            une_tol_f            % seuil de stationnarite des f_k             %
%            une_tol_g            % seuil validant EULER    en x_k             %
%                                   le sous-problème:
%                                   (1 exact, 2 pas de Cauchy, 3 CG tronqué)   %
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
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-SUPAERO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [                                                             ...
    x_opt  ,                                                           ...
    f_opt  ,                                                           ...
    fin    ,                                                           ...
    nit    ,                                                            ...
    x_seq]                                                                  ...
    = algo_SQP_BFGS(                                ...
    une_f          ,                                ...
    des_c          ,                                ...
    un_x0          ,                                ...
    un_lambda0     ,                                ...
    H0             ,                                ...
    une_tol_h      ,                                ...
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
gfdex           = (jac_fd(une_f,tempx)).';                                       
jcdex           = jac_fd(des_c,tempx)                                   ;
gldex = gfdex + jcdex'*templambda;
    
% Evaluer le Hessien de notre Lagrangian
H               = H0;
k               = 0;
fin             = 0;
Pzeros = zeros(p);

% first turn so everything is correctly initialized and we can compute
% difference
M = [H, jcdex'; jcdex, Pzeros];
B = [-gfdex; -cdex];
X = M\B; % d and lambda

x_seq = tempx';

while(fin==0)
    k = k + 1;
%     if isequal(k/10,round(k/10))
%         disp(k)
%     end
    dk = X(1:n); 
    templambda = X(n+1:end);
    tempx = tempx + dk;
    tempdx = dk;

    x_seq = [x_seq; tempx'];
    
    fdex = feval(une_f,tempx);
    cdex = feval(des_c,tempx);
    
    gfdex = (jac_fd(une_f,tempx)).';                                       
    jcdex = jac_fd(des_c,tempx);
   
    gldex_old = gldex;
    gldex = gfdex + jcdex'*templambda;
    ykm1 = gldex - gldex_old;
       
    if ykm1' * tempdx > 0
        H = H + ( ...
            ykm1 * ykm1')/(ykm1' * tempdx) - ( ...
            H * (tempdx * tempdx') * H)/(tempdx' * H * tempdx);
    end
    
    M = [H, jcdex'; jcdex, Pzeros];
    B = [-gfdex; -cdex];
    X = M\B; % d and lambda

    if (norm(gldex) < 0.01*une_tol_g) && (norm(dk) < 0.01*une_tol_x)
        fin = 2;
        % Insert patience criteria
    end    

    if(k==un_nit_max)
        fin =  3                                                           ;
    end
    
end
x_opt    =   tempx                                                         ;
f_opt    =    fdex                                                         ;

fprintf('The optimal value of x is X = (%6.4f, %6.4f) and the value of f is %6.4f.\n\n',x_opt(1),x_opt(2),f_opt)
nit      =   k                                                             ;

end

function jac = jac_fd(fun,X)
        X = X(:);
        N = length(X);
        S = zeros(N,9);
        funn = @(X)[ones(N,1).*fun(X)];
        
        if norm(X,2) < 0.1
            h = sqrt(eps);
        else
            h = sqrt(eps)*norm(X,2);
        end
        
        basis = eye(N).*h;
        C = [1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280];
       for j = 1:N
                for kk = 1:9
                    S(:,kk) = funn(X+(kk-5)*basis(:,j));
                end
       end
       jac = ((S*C')./h).';

end
function BE1_main(algo,prob)
addpath('./Problems');
warning off;
%**************************************************************************
%  Cette fonction constitue le noyau du programme principal relatif a
%  l'implementation de ce BE
%
%
%    algo: 1, 2, 3 ou 4 pour choisir l'algorithm.
%       algo=1 SQP
%       algo=2 SQP avec BFGS
%       algo=3 SQP avec BFGS et differences finies
%    prob: 1, 2 ou 3.
%
% Responsable: Y. Diouane (youssef.diouane@isae.fr) -- 
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-SUPAERO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%**************************
% GLOBALES EN MISE A JOUR *
%**************************

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

%***********************************************
%      Pour l'affichage (à ne pas modifier)     *
%***********************************************

ligne_tiret = ['-------------------------------------------------------' ]   ;
ligne_tiret = ['|' ligne_tiret  ligne_tiret  '|'];

lentete     = ['| DEPART|      METHODE      | PROBLEME |  FIN  |   F_COUNT   | ' ...
    ' G_COUNT   |  H_COUNT |   NITER    |   F_OPT    |'   ];

les_formats=['|%s| %s |%s |  %3.0f  |   %6.0f    |     %5.0f  |   %5.0f  |' ...
    '    %5.0f   |%10.5g  | '                      ];

nom_algo   ={ 'SQP avec derivees','  SQP avec BFGS  ' ,' SQP sans derivee'};
nom_prob   ={ '    P1   ','    P2   ' ,'   P3    '};

% en opti non-lineaire il faut un
% point de recherche initial ! x0
x(1,:)=[1 1];
x(2,:)=[2 2];
x(3,:)=[10 10];
x(4,:)=[2 5];
nom_point = {' (1,1) ',' (2,2) ','(10,10)',' (2,5) '};
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%      ETAPE DE RESOLUTION DU PROBLEME D'OPTIMISATION                     %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%
% Choix de probleme de minimisation (P1, P2 ou P3)
%
if(prob==1)
    f_fun    = @f_P1;
    g_fun    = @g_f_P1;
    h_fun    = @h_f_P1;
    c_fun    = @c_P1;
    jac_c_fun= @jac_c_P1;
    h_c_fun  = @h_c_P1;
    lambda0=[0];
    H0=eye(2);
elseif(prob==2)
    f_fun    = @f_P1;
    g_fun    = @g_f_P1;
    h_fun    = @h_f_P1;
    c_fun    = @c_P2;
    jac_c_fun= @jac_c_P2;
    h_c_fun  = @h_c_P2;
    lambda0=[0];
    H0=eye(2);
elseif(prob==3)
    f_fun=@f_rosenbrock;
    g_fun=@g_rosenbrock;
    h_fun=@h_rosenbrock;
    c_fun    = @c_P1;
    jac_c_fun= @jac_c_P1;
    h_c_fun  = @h_c_P1;
    lambda0=[0];
    H0=eye(2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Algorithme SQP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(algo==1)
    for i=1:4 % pour chaque point de part
        top_on(i)           =cputime;
        [x_opt(:,i),f_opt(i),conv(i),ite(i)]= ...
            algo_SQP( ...
            f_fun  ,...
            c_fun, ...
            x(i,:)          ,...
            lambda0, ...
            g_fun,...
            jac_c_fun, ...
            h_fun,...
            h_c_fun,...
            100,...
            10^-5,...
            10^-5,...
            10^-5 ...
            );
        
        % On garde les sorties pour l'afichage ..
        temps(i)                        =cputime - top_on(i)              ;
        x_optimale(:,i,algo)            =x_opt(:,i)                       ;
        t_cpu_time(algo ,i)             = temps(i)                        ;
        t_fin     (algo,i)              =conv(i)                          ;
        t_f_count (algo,i)              =f_count                          ;
        t_g_count (algo,i)              =g_count                          ;
        t_h_count (algo,i)              =h_count                          ;
        t_nit     (algo,i)              =ite(i)                           ;
        t_f_opt   (algo,i)              =f_opt(i)                         ;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Algorithme SQP avec BFGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(algo==2)
    for i=1:4 % pour chaque point de part
        top_on(i)           =cputime;
        [x_opt(:,i),f_opt(i),conv(i),ite(i)]= ...
            algo_SQP_BFGS( ...
            f_fun  ,...
            c_fun, ...
            x(i,:)          ,...
            lambda0, ...
            g_fun,...
            jac_c_fun, ...
            H0,...
            100,...
            10^-5,...
            10^-5,...
            10^-5 ...
            );
        
        % On garde les sorties pour l'afichage ..
        temps(i)                        =cputime - top_on(i)              ;
        x_optimale(:,i,algo)            =x_opt(:,i)                       ;
        t_cpu_time(algo ,i)             = temps(i)                        ;
        t_fin     (algo,i)              =conv(i)                          ;
        t_f_count (algo,i)              =f_count                          ;
        t_g_count (algo,i)              =g_count                          ;
        t_h_count (algo,i)              =h_count                          ;
        t_nit     (algo,i)              =ite(i)                           ;
        t_f_opt   (algo,i)              =f_opt(i)                         ;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Algorithme SQP sans derivee
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(algo==3)
    for i=1:4
        % pour chaque point de part
        top_on(i) =cputime;
        [x_opt(:,i),f_opt(i),conv(i),ite(i)]= ...
            algo_SQP_sans_derivee( ...
            f_fun,  ...
            c_fun,  ...
            x(i,:), ...
            lambda0,...
            H0     ,...
            10^-5  ,...
            100 ,...
            10^-5  ,...
            10^-5  ,...
            10^-5 ...
            );
        % On garde les sorties pour l'afichage ..
        temps(i)                        =cputime - top_on(i)              ;
        x_optimale(:,i,algo)            =x_opt(:,i)                       ;
        t_cpu_time(algo ,i)             = temps(i)                        ;
        t_fin     (algo,i)              =conv(i)                          ;
        t_f_count (algo,i)              =f_count                          ;
        t_g_count (algo,i)              =g_count                          ;
        t_h_count (algo,i)              =h_count                          ;
        t_nit     (algo,i)              =ite(i)                           ;
        t_f_opt   (algo,i)              =f_opt(i)                         ;
    end
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% ON AFFICHE TOUS LES RESULTATS  DE TOUS LES TESTS SOUS FORME D'UN TABLEAU%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
disp (sprintf     ('\n'                                                 ));
disp (ligne_tiret                                                        );
disp (lentete                                                            );
disp (ligne_tiret                                                        );

for i=1:4 % pour chaque point de départ
    x_optimale(:,i,algo);
    disp (sprintf     ( les_formats  ,                             ...
        nom_point    {i}    ,                                      ...
        nom_algo   {algo}   ,                                      ...
        nom_prob   {prob}   ,                                      ...
        t_fin      (algo,i) ,                                      ...
        t_f_count  (algo,i) ,                                      ...
        t_g_count  (algo,i) ,                                      ...
        t_h_count  (algo,i) ,                                      ...
        t_nit      (algo,i) ,                                      ...
        t_f_opt    (algo,i)  )            )                                  ;
end
disp (ligne_tiret                                                        );


clear           global      %on detruit les variables globales en sortant

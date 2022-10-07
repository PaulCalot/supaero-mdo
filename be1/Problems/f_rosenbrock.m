%*******************************************************************************
%   Codage de la fonction dite de "rosenbrock"                                 *
%*****************************************************************************0*
%                                                    %**************************
%                                                    % PARAMETRES EN ENTREE    *
%                                                    %**************************
%
%            x                    % un vecteur de R^2                          %
%
%                                                    %**************************
%                                                    % PARAMETRES EN SORTIE    *
%                                                    %**************************
%
%            fdex                 % rosenbrock(x)                              %
%
% Responsable: Y. Diouane (youssef.diouane@isae.fr) -- 2016/2017
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-Supaéro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function fdex = rosenbrock(x)        

%*****************************
%     LES DIFFERENTS "OBJETS"*
%*****************************
                                                     %**************************
                                                     % GLOBALES EN MISE A JOUR *
                                                     %**************************
          
          global f_count  ;                  % nombre     d'evaluations de 
                                             % fdex=rosenbrock(x) , sans pitie

                                                     %**************************
                                                     % FONCTIONS MATLAB        *
                                                     %**************************

%         return                             % utiliser help <nom_fonction>
%                      %                     % pour obtenir l'aide en ligne



%*********************
%      CODE          *
%*********************
 


fdex             = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2                          ;

f_count          = f_count   + 1                                               ;

return                                                                         ;

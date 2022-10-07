%*******************************************************************************
%   Codage du gradient de la fonction dite de "rosenbrock"                     *
%******************************************************************************0
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
%            gfdex                % rosenbrock'(x), un vecteur de R^2          %
%
% Responsable: Y. Diouane (youssef.diouane@isae.fr) -- 2016/2017
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-Supaéro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gfdex = g_rosenbrock(x)        

%*****************************
%     LES DIFFERENTS "OBJETS"*
%*****************************
                                                     %**************************
                                                     % GLOBALES EN MISE A JOUR *
                                                     %**************************

            global g_count  ;                % nombre     d'evaluations de
                                             % gfdex=rosenbrock'(x) , sans pitie
  
                                                     %**************************
                                                     % FONCTIONS MATLAB        *
                                                     %**************************

%           return                           % utiliser help <nom_fonction>
%                                            %   pour obtenir l'aide en ligne



%*********************
%      CODE          *
%*********************

gfdex            =[-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1))
 	            200*(x(2)-x(1)^2)                  ]                          ;
 
g_count          =       g_count          + 1                                  ;

return                                                                         ;

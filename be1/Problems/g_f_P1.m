%*******************************************************************************
%   Codage du gradient de la fonction objectif du probleme 1                   *
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
%            gfdex                % f_P1'(x), un vecteur de R^2          %
%
% Responsable: Y. Diouane (youssef.diouane@isae.fr) -- 2019/2020
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-SUPAERO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gfdex = g_f_P1(x)        

%*****************************
%     LES DIFFERENTS "OBJETS"*
%*****************************
                                                     %**************************
                                                     % GLOBALES EN MISE A JOUR *
                                                     %**************************

            global g_count  ;                % nombre     d'evaluations de
                                             % gfdex=bad'(x) , sans pitie
  
                                                     %**************************
                                                     % FONCTIONS MATLAB        *
                                                     %**************************

%           return                           % utiliser help <nom_fonction>
%                                            %   pour obtenir l'aide en ligne



%*********************
%      CODE          *
%*********************

gfdex            =2*[x(1)
 	               x(2)                      ]                                 ;
 
g_count          =       g_count          + 1                                  ;

return                                                                         ;

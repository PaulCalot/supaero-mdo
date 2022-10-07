%*******************************************************************************
%   Codage de la matrice Hessienne de la fonction objectif du probleme 1       *
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
%            hfdex                % h_f_P1"(x), une matrice de R^(2x2)     %
%
% Responsable: Y. Diouane (youssef.diouane@isae.fr) -- 2019/2020
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-SUAPERO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hfdex = h_f_P1(x)        

%*****************************
%     LES DIFFERENTS "OBJETS"*
%*****************************
                                                     %**************************
                                                     % GLOBALES EN MISE A JOUR *
                                                     %**************************

            global h_count  ;                % nombre     d'evaluations de
                                             % hfdex=h_f_P1"(x) , sans pitie
  
                                                     %**************************
                                                     % FONCTIONS MATLAB        *
                                                     %**************************

%           return                           % utiliser help <nom_fonction>
%                                            %   pour obtenir l'aide en ligne



%*********************
%      CODE          *
%*********************

hfdex=2*[1   0;0 1]                                          ;
h_count          =       h_count          + 1                                  ;

return                                                                         ;

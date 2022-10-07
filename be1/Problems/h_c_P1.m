%*******************************************************************************
%   Codage des Hessiens des contraintes du probleme 1       *
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
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-SUPAERO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hcdex = h_c_P1(x)        

%*****************************
%     LES DIFFERENTS "OBJETS"*
%*****************************
                                                     %**************************
                                                     % GLOBALES EN MISE A JOUR *
                                                     %**************************

            global hc_count  ;                % nombre     d'evaluations de
                                              % hcdex=c_P1''(x) , sans pitie
  
                                                     %**************************
                                                     % FONCTIONS MATLAB        *
                                                     %**************************

%           return                           % utiliser help <nom_fonction>
%                                            %   pour obtenir l'aide en ligne



%*********************
%      CODE          *
%*********************

hcdex=[0 0;0 0]                                                                ;
hc_count          =       hc_count          + 1                                ;

return                                                                         ;

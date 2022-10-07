%*******************************************************************************
%   Codage de la jacobien des contraintes du probleme 2                        *
%*******************************************************************************
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
%            jcdex                % c_P2'(x), une matrice de R^(1x2)     %
%
% Responsable: Y. Diouane (youssef.diouane@isae.fr) -- 2019/2020
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-Supaéro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function jcdex = jac_c_P2(x)        

%*****************************
%     LES DIFFERENTS "OBJETS"*
%*****************************
                                                     %**************************
                                                     % GLOBALES EN MISE A JOUR *
                                                     %**************************

            global jc_count  ;                % nombre     d'evaluations de
                                             % jcdex=c_P1'(x) , sans pitie
  
                                                     %**************************
                                                     % FONCTIONS MATLAB        *
                                                     %**************************

%           return                           % utiliser help <nom_fonction>
%                                            %   pour obtenir l'aide en ligne



%*********************
%      CODE          *
%*********************

jcdex=[2*x(1)*x(2) x(1)^2]                                                                    ;
jc_count          =       jc_count          + 1                                ;

return                                                                         ;

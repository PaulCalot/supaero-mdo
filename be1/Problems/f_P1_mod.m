%*******************************************************************************
%   Codage de la fonction objectif du probleme 1                               *
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
%            fdex                 % f_P1(x)                              %
%
% Responsable: Y. Diouane (youssef.diouane@isae.fr) -- 2019/2020
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-SUPAERO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function fdex = f_P1_mod(x)        

%*****************************
%     LES DIFFERENTS "OBJETS"*
%*****************************
                                                     %**************************
                                                     % GLOBALES EN MISE A JOUR *
                                                     %**************************
          
          global f_count  ;                  % nombre     d'evaluations de 
                                             % fdex=f_P1(x) , sans pitie

                                                     %**************************
                                                     % FONCTIONS MATLAB        *
                                                     %**************************

%         return                             % utiliser help <nom_fonction>
%                      %                     % pour obtenir l'aide en ligne



%*********************
%      CODE          *
%*********************
 


fdex             = x(:,1).^2 + x(:,2).^2                                             ;

f_count          = f_count   + 1                                               ;

return                                                                         ;

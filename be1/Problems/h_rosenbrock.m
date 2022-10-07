%*******************************************************************************
%   Codage de la matrice Hessienne de la fonction dite de "rosenbrock"         *
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
%            hfdex                % rosenbrock"(x), une matrice de R^(2x2)     %
%
% Responsable: Y. Diouane (youssef.diouane@isae.fr) -- 2016/2017
% (C) Institut Supérieur de l'Aéronautique et de l'Espace (ISAE-Supaéro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hfdex = h_rosenbrock(x)        

%*****************************
%     LES DIFFERENTS "OBJETS"*
%*****************************
                                                     %**************************
                                                     % GLOBALES EN MISE A JOUR *
                                                     %**************************

            global h_count  ;                % nombre     d'evaluations de
                                             % hfdex=rosenbrock"(x) , sans pitie
  
                                                     %**************************
                                                     % FONCTIONS MATLAB        *
                                                     %**************************

%           return                           % utiliser help <nom_fonction>
%                                            %   pour obtenir l'aide en ligne



%*********************
%      CODE          *
%*********************

hfdex=[1200 *x(1)^2-400*x(2)+2  -400*x(1) ;-400*x(1) 200]									;
h_count          =       h_count          + 1                                  ;

return                                                                         ;

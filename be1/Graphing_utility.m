%  Graphing utility

% First problem
% Things to plot
% contourplot of f
% series of point for every algorithm
% real minimum
clear xv yv
close all


clr = [0.9290 0.6940 0.1250];

        
% function_caller
prob = 1;
function_caller

    for i = 1:4
        figure(i)
        hold on
        arr = x_seq_SQP{i};
        a = max(abs(arr(:,1)));
        b = max(abs(arr(:,2)));

        xv = (-1.5*a:0.1:1.5*a)';
        yv = (-1.5*b:0.1:1.5*b)';
        c = 1 - xv; % Constraint
        [X, Y] = meshgrid(xv,yv);
        Z = f_P1_mod([X(:),Y(:)]);
        Z = reshape(Z,size(X));
        contour(X,Y,Z)
        
        plot(0.5,0.5,'Marker','pentagram','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','r')
        plot(x_seq_SQP{i}(1:end-1,1), x_seq_SQP{i}(1:end-1,2),'m-.o','MarkerFaceColor','m','MarkerSize',3,'LineWidth',1.2)
        plot(x_seq_SQP_BFGS{i}(1:end-1,1), x_seq_SQP_BFGS{i}(1:end-1,2),'g-.o','MarkerFaceColor','g','MarkerSize',3,'LineWidth',1.2)
        plot(x_seq_SQP_BFGS_SD{i}(1:end-1,1), x_seq_SQP_BFGS_SD{i}(1:end-1,2),'b-.o','MarkerFaceColor','b','MarkerSize',3,'LineWidth',1.2)
        plot(xv,c,'--','LineWidth',1.4,'Color',clr)

        plot(x_seq_SQP{i}(1,1), x_seq_SQP{i}(1,2),'Marker','diamond','MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','k')
        plot(x_seq_SQP{i}(end,1), x_seq_SQP{i}(end,2),'m-x','MarkerFaceColor','m','MarkerSize',10,'LineWidth',1.2)
        plot(x_seq_SQP_BFGS{i}(end,1), x_seq_SQP_BFGS{i}(end,2),'g-x','MarkerFaceColor','g','MarkerSize',10,'LineWidth',1.2)
        plot(x_seq_SQP_BFGS_SD{i}(end,1), x_seq_SQP_BFGS_SD{i}(end,2),'b-x','MarkerFaceColor','b','MarkerSize',10,'LineWidth',1.2)
        
        legl = {'f(X)','Target','SQP','BFGS',' BFGS - FD','c(X)','Start'};
        leg = legend(legl,'Location','bestoutside');

        axis([-1.5*a 1.5*a -1.5*b 1.5*b])
        xlabel_1 = xlabel('x');
        ylabel_1 = ylabel('y');
        tl = ['Problem $P_1$ with $ \mathbf{X_0} = [',num2str(x_seq_SQP{i}(1,1)),',',num2str(x_seq_SQP{i}(1,2)),']$'];
        ttl = title(tl);
        grid minor
        
        set(leg,   'Interpreter', 'latex', 'fontsize', 12);
        set(xlabel_1,'Interpreter', 'latex', 'fontsize', 15, 'LineWidth', 1.5 );
        set(ylabel_1,'Interpreter', 'latex', 'fontsize', 15, 'LineWidth', 1.5 );
        set(ttl, 'Interpreter', 'latex', 'fontsize', 16, 'LineWidth', 1.5);
    end

clear xv yv
prob = 2;
function_caller    
    
    for i = 1:4
        figure(4+i)
        hold on
        arr = x_seq_SQP{i};
        a = max(abs(arr(:,1)));
        b = max(abs(arr(:,2)));

        xv = (-1.5*a:0.1:1.5*a)';
        yv = (-1.5*b:0.1:1.5*b)';
        c = 16./xv.^2; % Constraint
        [X, Y] = meshgrid(xv,yv);
        Z = f_P1_mod([X(:),Y(:)]);
        Z = reshape(Z,size(X));
        contour(X,Y,Z)
        plot(2.52,2.52,'Marker','pentagram','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','r')

        plot(x_seq_SQP{i}(1:end-1,1), x_seq_SQP{i}(1:end-1,2),'mo-.','MarkerFaceColor','m','MarkerSize',3,'LineWidth',1.2)
        plot(x_seq_SQP_BFGS{i}(1:end-1,1), x_seq_SQP_BFGS{i}(1:end-1,2),'go-.','MarkerFaceColor','g','MarkerSize',4,'LineWidth',1.2)
        plot(x_seq_SQP_BFGS_SD{i}(1:end-1,1), x_seq_SQP_BFGS_SD{i}(1:end-1,2),'bo-.','MarkerFaceColor','b','MarkerSize',5,'LineWidth',1.2)
        plot(xv,c,'--','LineWidth',1.4,'Color',clr)
        
        plot(x_seq_SQP{i}(1,1), x_seq_SQP{i}(1,2),'Marker','diamond','MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','k')
        plot(x_seq_SQP{i}(end,1), x_seq_SQP{i}(end,2),'m-x','MarkerFaceColor','m','MarkerSize',10,'LineWidth',1.2)
        plot(x_seq_SQP_BFGS{i}(end,1), x_seq_SQP_BFGS{i}(end,2),'g-x','MarkerFaceColor','g','MarkerSize',10,'LineWidth',1.2)
        plot(x_seq_SQP_BFGS_SD{i}(end,1), x_seq_SQP_BFGS_SD{i}(end,2),'b-x','MarkerFaceColor','b','MarkerSize',10,'LineWidth',1.2)
        
        legl = {'f(X)','Target','SQP','BFGS',' BFGS - FD','c(X)','Start'};
        leg = legend(legl,'Location','bestoutside');

        axis([-1.5*a 1.5*a -1.5*b 1.5*b])
        xlabel_1 = xlabel('x');
        ylabel_1 = ylabel('y');
        tl = ['Problem $P_2$ with $ \mathbf{X_0} = [',num2str(x_seq_SQP{i}(1,1)),',',num2str(x_seq_SQP{i}(1,2)),']$'];
        ttl = title(tl);
        grid minor
        
        set(leg,   'Interpreter', 'latex', 'fontsize', 12);
        set(xlabel_1,'Interpreter', 'latex', 'fontsize', 15, 'LineWidth', 1.5 );
        set(ylabel_1,'Interpreter', 'latex', 'fontsize', 15, 'LineWidth', 1.5 );
        set(ttl, 'Interpreter', 'latex', 'fontsize', 16, 'LineWidth', 1.5);
    end    

prob = 3;
function_caller 
clear xv yv
    
    for i = 1:4
        figure(8+i)
        hold on
        arr = x_seq_SQP{i};
        a = max(abs(arr(:,1)));
        b = max(abs(arr(:,2)));

        xv = (-1.5*a:0.1:1.5*a)';
        yv = (-1.5*b:0.1:1.5*b)';
        c = 1 - xv; % Constraint
        [X, Y] = meshgrid(xv,yv);
        Z = rosenbrock_mod([X(:),Y(:)]);
        Z = reshape(Z,size(X));
        contour(X,Y,Z)
        
        plot(0.5,0.5,'Marker','pentagram','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','r')        
        plot(x_seq_SQP{i}(1:end-1,1), x_seq_SQP{i}(1:end-1,2),'mo-.','MarkerFaceColor','m','MarkerSize',3,'LineWidth',1.2)
        plot(x_seq_SQP_BFGS{i}(1:end-1,1), x_seq_SQP_BFGS{i}(1:end-1,2),'go-.','MarkerFaceColor','g','MarkerSize',4,'LineWidth',1.2)
        plot(x_seq_SQP_BFGS_SD{i}(1:end-1,1), x_seq_SQP_BFGS_SD{i}(1:end-1,2),'bo-.','MarkerFaceColor','b','MarkerSize',5,'LineWidth',1.2)
        plot(xv,c,'--','LineWidth',1.4,'Color',clr)

        
        plot(x_seq_SQP{i}(1,1), x_seq_SQP{i}(1,2),'Marker','diamond','MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','k')
        plot(x_seq_SQP{i}(end,1), x_seq_SQP{i}(end,2),'m-x','MarkerFaceColor','m','MarkerSize',10,'LineWidth',1.2)
        plot(x_seq_SQP_BFGS{i}(end,1), x_seq_SQP_BFGS{i}(end,2),'g-x','MarkerFaceColor','g','MarkerSize',10,'LineWidth',1.2)
        plot(x_seq_SQP_BFGS_SD{i}(end,1), x_seq_SQP_BFGS_SD{i}(end,2),'b-x','MarkerFaceColor','b','MarkerSize',10,'LineWidth',1.2)        
        legl = {'f(X)','Target','SQP','BFGS',' BFGS - FD','c(X)','Start'};
        leg = legend(legl,'Location','bestoutside');

        axis([-1.5*a 1.5*a -1.5*b 1.5*b])
        xlabel_1 = xlabel('x');
        ylabel_1 = ylabel('y');
        tl = ['Problem $P_3$ with $ \mathbf{X_0} = [',num2str(x_seq_SQP{i}(1,1)),',',num2str(x_seq_SQP{i}(1,2)),']$'];
        ttl = title(tl);
        grid minor

        set(leg,   'Interpreter', 'latex', 'fontsize', 12);
        set(xlabel_1,'Interpreter', 'latex', 'fontsize', 15, 'LineWidth', 1.5 );
        set(ylabel_1,'Interpreter', 'latex', 'fontsize', 15, 'LineWidth', 1.5 );
        set(ttl, 'Interpreter', 'latex', 'fontsize', 16, 'LineWidth', 1.5);
    end    

    
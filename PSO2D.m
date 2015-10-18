
% SIMULATING HYBRID BOIDS-PSO ALGORITHM FOR MULTI-ROBOT SYSTEMS
% IN UNKNOWN ENVIRONMENT EXPLORATION
%
% THESIS IN ELECTRONICS AND TELECOMMUNICATIONS MAJOR
% AUTHOR: HOANG ANH QUY
% SUPERVISOR: Dr. PHAM MINH TRIEN
%
% FACULTY OF ELECTRONICS AND TELECOMMUNICATIONS
% UNIVERSITY OF ENGINEERING AND TECHNOLOGY
% VIETNAM NATIONAL UNIVERSITY, HANOI
%
% May 2015



%//////////////////////////////////////////////////////////////////////////
% INITIALIZATION
%//////////////////////////////////////////////////////////////////////////


clc;
close all; clear all; 

set(0, 'DefaultFigurePosition', [299, 0, 768, 768]);
color = [1 0 0];
particleColor = [0 0 1];
WAIT = .2;



global noParticle range spaceSize vLimitX vLimitY maxStep commuRange
noParticle = 15; % number of particles
maxStep = 200; % number of iterations
range = 1; % physical dimension of a particle
sensingRange = range*ones(1,noParticle);
convergence = 10000;

% staticObs = rand(2,3);
staticObs = [-35 -20; 35 -30];
 
movingObs = [5;-35]; direction = 1;


% Coefficients ------------------------------------------------------------
w = linspace(1,1,maxStep); % Inertial factor
c1 = linspace(2,0.1,maxStep); % Personal impact factor
c2 = linspace(0.1,2,maxStep); % Social impact factor

% Space initialization ----------------------------------------------------
upLimit_1 = 50; dwLimit_1 = -50; assert(dwLimit_1 < upLimit_1); %1: x, 2: y
upLimit_2 = 50; dwLimit_2 = -50; assert(dwLimit_2 < upLimit_2);

% Virtual limit for avoiding collision with boundaries --------------------
dwSpace = repmat([dwLimit_1; dwLimit_2],1,noParticle);
upSpace = repmat([upLimit_1; upLimit_2],1,noParticle);

% Virtual limits for avoiding collision with boundaries, regarding size of particles
dw = dwSpace + range;
up = upSpace - range; % the swarm's coverage is almost all the space, except for four corners.
spaceSize = upSpace - dwSpace;

swarm = dw + rand(2,noParticle).*spaceSize; % initial position ------------
% swarm = initialSwarm - 40; % evenly distributed after run BOIDS-PSOA
% swarm = matranGoc;

veloScale = 0.5; % the ratio of maximum velocity over particle size -------
vLimitX = veloScale*range;
vLimitY = veloScale*range;
velo = repmat([vLimitX; vLimitY],1,noParticle); % initial velocity

commuRange = 10*range; % communication range, connected distance ----------

%//////////////////////////////////////////////////////////////////////////










% % CONNECTIVITY CHECK-------------------------------------------------------
%
% % Adjacency Matrix
% A = zeros(noParticle,noParticle);
% for i1 = 1:noParticle
%     for i2 = (i1+1):noParticle
%         if r(swarm(:,i1),swarm(:,i2)) < commuRange
%             A(i1,i2) = 1;
%             A(i2,i1) = 1;
%         end
%     end
%     A(i1,i1) = 0;
% end
%
% % Degree Matrix
% D = zeros(noParticle,noParticle);
% for i1 = 1:noParticle
%     for i2 = 1:noParticle
%         if (r(swarm(:,i1),swarm(:,i2)) < commuRange) && (i1 ~= i2)
%             D(i1,i1) = D(i1,i1)+1;
%         end
%     end
% end
%
% % Laplacian Matrix
% L = D - A;
%
% eigen = sort(eig(L));
% lambda2 = zeros(1,maxStep);
%
% % -------------------------------------------------------------------------










%//////////////////////////////////////////////////////////////////////////
% OBJECTIVE FUNCTION
%//////////////////////////////////////////////////////////////////////////


% PLOT CONTOUR OF OBJECTIVE FUNCTION --------------------------------------

% -----------------THREE-HUMP CAMEL FUNCTION 3D PLOT -----------------------
% x = [-2:.1:2];
x = [-50:.1:50];
y = x;
[xx,yy] = meshgrid(x,y);
% zz = 2*xx.^2 - 1.05*xx.^4 + xx.^6/6 + yy.^2 + yy.*xx;
zz = 2*(xx/25).^2 - 1.05*(xx/25).^4 + (xx/25).^6/6 + (yy/25).^2 + (yy/25).*(xx/25);
% -----------------------------------------------------------------------

% %-----------------BOHACHEVSKY FUNCTION 3D PLOT --------------------------
% x = [-50:.1:50];
% y = x;
% [xx,yy] = meshgrid(x,y);
% zz = xx.^2 + 2*yy.^2 - 0.3*cos(3*pi*xx + 4*pi*yy) + 0.3;
% % -----------------------------------------------------------------------

% %-----------------SPHERE FUNCTION 3D PLOT -------------------------------
% x = [-50:.1:50];
% y = x;
% [xx,yy] = meshgrid(x,y);
% zz = xx.^2 + yy.^2;
% % -----------------------------------------------------------------------

% %-----------------ROSENBROCK FUNCTION 3D PLOT -----------------------------
% x = [-50:.1:50];
% y = x;
% k = 25;
% [xx,yy] = meshgrid(x,y);
% zz = ((xx+k)/k-1).^2 + 100*((yy+k)/k-((xx+k)/k).^2).^2;
% % -----------------------------------------------------------------------

% -------------------------------------------------------------------------










% MANIPULATE FITNESSES ----------------------------------------------------
fit = subFitness(swarm);

swarmX = swarm(1,:);
swarmY = swarm(2,:);

[gBestVal,gBestIndex] = min(fit);
gBest = repmat(swarm(:,gBestIndex),1,noParticle);

pBest = swarm; % initial personal best
pBestFit = subFitness(pBest);

BestFitnessEver = zeros(1,maxStep);
MeanFitnessEver = zeros(1,maxStep);
%--------------------------------------------------------------------------










%//////////////////////////////////////////////////////////////////////////
% UPDATING and TERMINATING SEARCHING PROCESS
%//////////////////////////////////////////////////////////////////////////

for step = 1:maxStep
    % UPDATE VELOCITY /////////////////////////////////////////////////////
    
%         % Conventional PSO velocity
%         velo = (w(step)*velo + c1(step)*rand*(pBest-swarm) + c2(step)*rand*(gBest-swarm));
%     
    % BOIDS-PSOA velocity
    velo = (w(step)*velo + c1(step)*rand*(pBest-swarm) + c2(step)*rand*(gBest-swarm))...
        + updateBoidVelocity(swarm,velo,staticObs,movingObs);
%     
    
    % Check speed limit ---------------------------------------------------
    for parNo = 1:noParticle % no of particle
        s = velo(1,parNo)/velo(2,parNo); % ratio of v components
        sign1 =  velo(1,parNo)/abs(velo(1,parNo));
        sign2 =  velo(2,parNo)/abs(velo(2,parNo));
        if (abs(velo(1,parNo)) > vLimitX)
            velo(1,parNo) = rand*vLimitX*sign1;
            velo(2,parNo) = velo(1,parNo)/s;
        elseif (abs(velo(2,parNo)) > vLimitY)
            velo(2,parNo) = rand*vLimitY*sign2;
            velo(1,parNo) = velo(2,parNo)*s;
        end
    end
    % ---------------------------------------------------------------------
    
    
    
    
    % UPDATE POSITION /////////////////////////////////////////////////////
    tempPos = swarm + velo;
    % Check space limit
    for iDim = 1:2 %dimension index
        for parNo = 1:noParticle %no of particle
            if tempPos(iDim,parNo) > up(iDim,parNo)
                tempPos(iDim,parNo) = up(iDim,parNo) - rand*(up(iDim,parNo) - swarm(iDim,parNo));
            elseif tempPos(iDim,parNo) < dw(iDim,parNo)
                tempPos(iDim,parNo) = dw(iDim,parNo) - rand*(dw(iDim,parNo) - swarm(iDim,parNo));
            end
        end
    end
    swarm = tempPos;
    % /////////////////////////////////////////////////////////////////////
    
    
    % UPDATE POSITION OF MOVING OBSTACLES
    if movingObs(2) >= 35
        direction = -1;
    elseif movingObs(2) <= -35
        direction = 1;
    end
    
    if direction == 1
        movingObs(2) = movingObs(2) + vLimitY/2;
    else 
        movingObs(2) = movingObs(2) - vLimitY/2;
    end
     
    
    
    
    
    
    
    
    
    
    % UPDATE BEST POSITION AND LEADER /////////////////////////////////////
    
    % Find new personal best pBestFit = subFitness(pBest);
    newFit = fitness(swarm); % sensing range is applied
    % newFit = subFitness(swarm); % without sensing range
    pBestIndex = find(newFit<pBestFit);
    pBestFit(:,pBestIndex) = newFit(:,pBestIndex); % update best personal value
    pBest(:,pBestIndex) = swarm(:,pBestIndex); % update best personal position
    
    % Find new global best
    [gBestValTemp,gBestIndexTemp] = min(newFit);
    if gBestValTemp < gBestVal
        gBestVal = gBestValTemp;
        gBestIndex = gBestIndexTemp;
        gBest = repmat(swarm(:,gBestIndex),1,noParticle);
    end
    BestFitnessEver(step) = gBestVal;
    MeanFitnessEver(step) = mean(newFit);
    
    fit = [fit;newFit]; % fitness matrix
    swarmX = [swarmX; swarm(1,:)]; % x-coordinator (of all steps) matrix
    swarmY = [swarmY; swarm(2,:)]; % y-coordinator matrix
    %----------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    
%     % CONNECTIVITY CHECK //////////////////////////////////////////////////
%     % Adjacency Matrix
%     A = zeros(noParticle,noParticle);
%     for i1 = 1:noParticle
%         for i2 = (i1+1):noParticle
%             if r(swarm(:,i1),swarm(:,i2)) < commuRange
%                 A(i1,i2) = 1;
%                 A(i2,i1) = 1;
%             end
%         end
%         A(i1,i1) = 0;
%     end
%     
%     % Degree Matrix
%     D = zeros(noParticle,noParticle);
%     for i1 = 1:noParticle
%         for i2 = 1:noParticle
%             if (r(swarm(:,i1),swarm(:,i2)) < commuRange) && i1 ~= i2
%                 D(i1,i1) = D(i1,i1)+1;
%             end
%         end
%         %      D(i1,i1) =  D(i1,i1) - 1; % remove the increase in degree by itself
%     end
%     
%     % Laplacian Matrix
%     L = D - A;
%     
%     eigen = sort(eig(L));
%     lambda2(step) = eigen(2);
%    
%     % -------------------------------------------------------------------------   
    
    
    
    
    
    
    
    
    % CHECK TERMINATING CONDITION /////////////////////////////////////////
    %     distance = 0;
    %     for k = 2:noParticle
    %         ds = sqrt((swarm(1,1)-swarm(1,k))^2 + (swarm(2,1)-swarm(2,k))^2 );
    %         if(ds > distance)
    %             distance = ds;
    %         end
    %     end
    %     if (pi*distance^2 < 0.001*(upLimit_1 - dwLimit_1)*(upLimit_2 - dwLimit_2))
    %         break;
    %     end
    %     convergence = distance;
    
    
    
    
    
    
    
    
    
    
    % GRAPHIC /////////////////////////////////////////////////////////////
    
    % Motion tracking -----------------------------------------------------
    clf; % clear current figure
%     figure('position', [0, 0, 768, 768]);

    scatter(0,0,100,color,'filled');
%     scatter(0,0,'rO','filled'); % highlight target
    hold on

    scatter(staticObs(1,:),staticObs(2,:),100,[0 0 0],'filled'); % draw static obstacles
    hold on
    scatter(movingObs(1),movingObs(2),100,[0 1 0.5],'filled');
    hold on
    
    
    xlabel('x coordinator','FontSize',16)
    ylabel('y coordinator','FontSize',16)
    a = [num2str(step)];
    str = ['Step: ',num2str(step),'\it \color[rgb]{0.2,0.5,0.5} Best Value = ', num2str(gBestVal)];
    title(str,'FontSize',18);
    set(gca,'fontsize',16);
    hold on
    
%     scatter(swarmX(step,:),swarmY(step,:),'r+') % draw all the swarm
%     hold on
%     
%     scatter(swarmX(step,:),swarmY(step,:),500,color);
%     
    swarmStep = [swarmX(step,:);swarmY(step,:)]';
%     viscircles(swarmStep,sensingRange,'EdgeColor','b'); 
    scatter(swarmX(step,:),swarmY(step,:),90,particleColor,'filled'); % visual alternation of the above line
    % draw all the swarm with each particle's size
    
%     viscircles(swarmStep,commuRange*sensingRange,'EdgeColor','g','LineWidth',0.1);
%     % draw communication range

 viscircles(swarmStep,8*sensingRange,'EdgeColor','g','LineWidth',0.1);
%     % draw collision range
    
    axis([dwLimit_1-1 upLimit_1+1 dwLimit_2-1 upLimit_2+1]);
    hold on;
    
    
    %     str2 = {['\lambda_2 = ',num2str(lambda2(step))]};
    %     annotation('textbox', [0.66,0.7,0.1,0.1],'String', str2) % add notation on graph
    
    pause(WAIT/5);
    
    % Get images for report
%     if step == 1 || step == 5 || step == 20 || step == 40 || step == 60 || step == 100 || step == 150 || step == 200
    if step == 1 || step == 200  
        contour(xx,yy,zz,40);
%         saveas(gcf,num2str(step),'emf'); % save current figure in a defined format
    end
    
    %
    %     % Trajectory tracking
    %     ------------------------------------------------- scatter(0,0,'bo');
    %     xlabel('x coordinator','FontSize',16) ylabel('y
    %     coordinator','FontSize',16) a = [num2str(step)]; str = ['Step:
    %     ',num2str(step),'\it \color[rgb]{0.2,0.5,0.5} Best Value = ',
    %     num2str(gBestVal)]; title(str,'FontSize',12) hold on
    %     scatter(swarmX(step,:),swarmY(step,:)) axis([dwLimit_1-1 upLimit_1+1
    %     dwLimit_2-1 upLimit_2+1]); hold on; pause(WAIT/5);
    
    
end

% initialSwarm = swarm; % create a uniformly distributed swarm










%//////////////////////////////////////////////////////////////////////////
% DISPLAY RESULTS
%//////////////////////////////////////////////////////////////////////////

pause(5*WAIT)
hold off

figure;

% subplot(3,1,1); plot(BestFitnessEver,'color','r','LineWidth',2);
% title311 = ['Best fitness (final value = ',num2str(gBestVal),' )'];
% grid;title(['\fontsize{12} {\color{blue}Best fitness (final value = ',num2str(gBestVal),' )}']);
% xlabel('Step ','FontSize',10)
% ylabel('Fitness','FontSize',10)
% 
% subplot(3,1,2); plot(MeanFitnessEver,'LineWidth',2);
% grid;title(['\fontsize{12} {\color{blue}Mean of fitness}']);
% xlabel('Step ','FontSize',10)
% ylabel('Fitness','FontSize',10)
% 
% subplot(3,1,3);

% plot(BestFitnessEver,'color','r','LineWidth',3);
% hold on; plot(MeanFitnessEver,'LineWidth',3);
% grid; title(['\fontsize{18} {\color{blue}Mean of fitness \color{red}vs \color{blue}Best fitness}']);
% xlabel('Step ','FontSize',16)
% ylabel('Fitness','FontSize',16)
% legend('Best Fitness','MeanFitness');
% set(gca,'fontsize',16);
% 
% saveas(gcf,'Final result','emf');

pause(10*WAIT)

% figure;
% plot(lambda2,'LineWidth',2);
% grid;title('Algebraic Connectivity','Fontsize',20);
% xlabel('Step','FontSize',15)
% ylabel('\lambda_2','FontSize',15)
% saveas(gcf,'Algebraic Connectivity','emf');
% % axis([0 maxStep -1 4]);

result = ['Best value is ',num2str(gBestVal),', found at (',num2str(gBest(1,1)),',',num2str(gBest(2,1)),')'];
disp(result);

% END //////////////////////////////////////////////////////////////////////










% MFILE TO DRAW OBJECTIVE FUNCTIONS ///////////////////////////////////////

% x = [-50:.1:50];
% y = x;
% [xx,yy] = meshgrid(x,y);
% zz = xx.^2 + 2*yy.^2 - 0.3*cos(3*pi*xx + 4*pi*yy) + 0.3;

% %-----------------BOHACHEVSKY FUNCTION 3D PLOT --------------------------
% x = [-50:.1:50];
% y = x;
% [xx,yy] = meshgrid(x,y);
% zz = xx.^2 + 2*yy.^2 - 0.3*cos(3*pi*xx + 4*pi*yy) + 0.3;
% % -----------------------------------------------------------------------

% %-----------------SPHERE FUNCTION 3D PLOT -------------------------------
% x = [-50:.1:50];
% y = x;
% [xx,yy] = meshgrid(x,y);
% zz = xx.^2 + yy.^2;
% % -----------------------------------------------------------------------

%-----------------ROSENBROCK FUNCTION 3D PLOT -------------------------------
% x = [-50:.1:50];
% y = x;
% k = 25;
% [xx,yy] = meshgrid(x,y);
% zz = ((xx+k)/k-1).^2 + 100*((yy+k)/k-((xx+k)/k).^2).^2;
% -----------------------------------------------------------------------

% % -----------------THREE-HUMP CAMEL FUNCTION 3D PLOT -----------------------
% % x = [-2:.1:2];
% k = 25;
% x = [-50:.1:50];
% y = x;
% [xx,yy] = meshgrid(x,y);
% zz = 2*(xx/k).^2 - 1.05*(xx/k).^4 + (xx/k).^6/6 + (yy/k).^2 + (yy/k).*(xx/k);
% % -------------------------------------------------------------------------

% contour(xx,yy,zz,40);
% xlabel('x coordinator ','FontSize',10)
% ylabel('y coordinator','FontSize',10)

% surf(xx,yy,zz)
% axis tight
% shading interp
% colorbar
% xlabel('x coordinator','FontSize',10)
% ylabel('y coordinator','FontSize',10)
% zlabel('Fitness','FontSize',10)
% 
% saveas(gcf,'Rosen compare surf','emf');
% close all;

%//////////////////////////////////////////////////////////////////////////









% MFILE TO CALCULATE FITNESS FUNCTION /////////////////////////////////////

% function x = fitness(swarm)
% 
% global noParticle
% global range
% 
% diagonal = range/sqrt(2);
% rg = range*ones(1,noParticle);
% dia = diagonal*ones(1,noParticle);
% zr = zeros(1,noParticle);
% 
% % swarm = dwSpace + rand*spaceSize; % initial position - real position 
% swarmN  = swarm + [zr; rg]; % x, y
% swarmNW = swarm + [dia; dia];
% swarmW  = swarm + [rg; zr];
% swarmSW = swarm + [dia; -dia];
% swarmS  = swarm + [zr; -rg];
% swarmSE = swarm + [-dia; -dia];
% swarmE  = swarm + [-rg; zr];
% swarmNE = swarm + [-dia; dia];
% 
% subN = subFitness(swarmN);
% subNW = subFitness(swarmNW);
% subW = subFitness(swarmW);
% subSW = subFitness(swarmSW);
% subS = subFitness(swarmS);
% subSE = subFitness(swarmSE);
% subE = subFitness(swarmE);
% subNE = subFitness(swarmNE);
% sub = subFitness(swarm);
% 
% tempFitness = [sub; subN; subNW; subW; subSW; subS; subSE; subE; subNE];
% x = min(tempFitness);
% 
% %//////////////////////////////////////////////////////////////////////////










% MFILE TO CALCULATE SUBFITNESS
% 
% function sub = subFitness(swarm)
% global noParticle
% sub = ones(1,noParticle);

% %-----------------BOHACHEVSKY FUNCTION
% for i = 1:noParticle
%     sub(i) = swarm(1,i)^2 + 2*swarm(2,i)^2 - 0.3*cos(3*pi*swarm(1,i) + 4*pi*swarm(2,i)) + 0.3;
% end
%--------------------------------------------------------------------------


% %-----------------THREE-HUMP CAMEL FUNCTION -----------------------------
% k = 25;
% for i = 1:noParticle
% 	sub(i) = 2*(swarm(1,i)/k)^2 - 1.05*(swarm(1,i)/k)^4 + (swarm(1,i)/k)^6/6 + (swarm(1,i)/k)*(swarm(2,i)/k) + (swarm(2,i)/k)^2;
% end
%--------------------------------------------------------------------------


% %-----------------SPHERE FUNCTION----------------------------------------
% sub = (sum(swarm.*swarm));
%--------------------------------------------------------------------------


% %-----------------ROSENBROCK FUNCTION -------------------------------------
% k = 25;
% for i = 1:noParticle
%    sub(i) = ((swarm(1,i)+k)/k-1).^2 + 100*((swarm(2,i)+k)/k-((swarm(1,i)+k)/k).^2).^2;
% end
% %--------------------------------------------------------------------------

%//////////////////////////////////////////////////////////////////////////










% FUNCTION TO IMPLEMENT BOIDS
%
% function veloOutput = updateBoidVelocity(swarm,velo,convergence)
% 
% global noParticle
% global range
% 
% veloOutput = zeros(2,noParticle);
% mirrorVelo = velo;
% vSeparation = zeros(2,noParticle);
% vAlignment = zeros(2,noParticle);
% vCohesion = zeros(2,noParticle);
% centerCohe = zeros(2,noParticle); % notice this point
% neighbor = 0;
% 
% rSeparation = 8*range;
% % rSeparation = range + convergence/3;
% rCohesion = 20*range;
% 
% % Update vSeparation
% for boid1 = 1:noParticle
%     for boid2 = 1:noParticle
%         if boid2 ~= boid1
%             if r(swarm(:,boid1),swarm(:,boid2)) < rSeparation
%                 r12 = r(swarm(:,boid1),swarm(:,boid2));
%                 magnitude = rSeparation - r12;
%                 connect = swarm(:,boid1) - swarm(:,boid2);
%                 vX = magnitude/r12*connect(1,1);
%                 vY = magnitude/r12*connect(2,1);
%                 vSeparation(:,boid1) = vSeparation(:,boid1) + [vX;vY];
%             end
%         end
%     end
% end
% 
% % Update vCohesion
% for boid1 = 1:noParticle
%     for boid2 = 1:noParticle
%         if boid2 ~= boid1
%             if r(swarm(:,boid1),swarm(:,boid2)) < rCohesion
%                 neighbor = neighbor + 1;
%                 centerCohe(:,boid1) = centerCohe(:,boid1) + swarm(:,boid2);
%             end
%         end
%     end
%     if neighbor ~= 0
%         centerCohe(:,boid1) = centerCohe(:,boid1)/neighbor;
%         neighbor = 0;
%         vCohesion(:,boid1) = (centerCohe(:,boid1) - swarm(:,boid1))/100;
%         % move 1% towards the center of neighbors
%     end
% end
% 
% neighbor = 0;
% % Update vAlignment
% for boid1 = 1:noParticle
%     for boid2 = 1:noParticle
%         if boid2 ~= boid1
%             if r(swarm(:,boid1),swarm(:,boid2)) < 20*range
%                 mirrorVelo(:,boid1) = mirrorVelo(:,boid1) + velo(:,boid2);
%                 neighbor = neighbor + 1;
%             end
%         end
%     end
%     if neighbor ~= 0
%         vAlignment(:,boid1) = mirrorVelo(:,boid1)/neighbor;
%         neighbor = 0;
%         vAlignment(:,boid1) = (vAlignment(:,boid1) - velo(:,boid1))/10;
%     end
% end
% 
% veloOutput = velo + vSeparation*2 + vCohesion*6 + vAlignment*2;
% end
% 
% % function sqr = sqr(x)
% % sqr = x.*x;
% % end
% % 
% % function distance = r(a,b)
% % distance = sqrt(sum(sqr(a-b)));
% % end

%//////////////////////////////////////////////////////////////////////////












function HBPSO(objectiveFunction,connectivityCheck)

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
% COMPLETED IN MAY 2015
% 
% objectiveFunction is one of the following: threeHumpCamel (Three-Hump Camel Function),
% bohachevsky (Bohachevsky Function), sphere (Sphere Function), rosenbrock (Rosenbrock Function)
% connectivityCheck is either 1 (show calculation for lambda2) or 0




% |INITIALIZATION||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
clc; close all; %clear all; % clear previous data, clear screen
set(0, 'DefaultFigurePosition', [299, 0, 768, 768]); % make every graph fit my screen (1366*768)
red = [1 0 0]; green = [0 1 0]; blue = [0 0 1]; black = [0 0 0]; % set colour for later use
WAIT = .2; % interval to display steps

global noParticle range spaceSize vLimitX vLimitY maxStep commuRange
noParticle = 15; % number of particles
maxStep = 200; % number of iterations
range = 1; % physical radius of a particle
sensingRange = range*ones(1,noParticle);

% staticObs = rand(2,3);
staticObs = [-35 -21; 35 -31]; % Position of static obstacles
movingObs = [5;-35]; direction = 1; % Position and moving direction of moving obstacles


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

if connectivityCheck == 0
    swarm = dw + rand(2,noParticle).*spaceSize; % initial position --------
elseif connectivityCheck == 1
    load('arrangePost.mat'); swarm = swarm - 30;
else disp('Invalid connection property!');
end

veloScale = 0.5; % the ratio of maximum velocity over particle size -------
vLimitX = veloScale*range;
vLimitY = veloScale*range;
velo = repmat([vLimitX; vLimitY],1,noParticle); % initial velocity

commuRange = 10*range; % communication range, connected distance ----------
lambda2 = zeros(1,maxStep);




% | OBJECTIVE FUNCTION ||||||||||||||||||||||||||||||||||||||||||||||||||||||

% PLOT CONTOUR OF OBJECTIVE FUNCTION --------------------------------------
switch objectiveFunction
    case 'threeHumpCamel'
        % -----------------THREE-HUMP CAMEL FUNCTION 3D PLOT --------------
        % x = [-2:.1:2];
        x = [-50:.1:50];
        y = x;
        [xx,yy] = meshgrid(x,y);
        % zz = 2*xx.^2 - 1.05*xx.^4 + xx.^6/6 + yy.^2 + yy.*xx;
        zz = 2*(xx/25).^2 - 1.05*(xx/25).^4 + (xx/25).^6/6 + (yy/25).^2 + (yy/25).*(xx/25);
        
    case 'bohachevsky'
        %-----------------BOHACHEVSKY FUNCTION 3D PLOT --------------------
        x = [-50:.1:50];
        y = x;
        [xx,yy] = meshgrid(x,y);
        zz = xx.^2 + 2*yy.^2 - 0.3*cos(3*pi*xx + 4*pi*yy) + 0.3;
        
    case 'sphere'
        %-----------------SPHERE FUNCTION 3D PLOT -------------------------
        x = [-50:.1:50];
        y = x;
        [xx,yy] = meshgrid(x,y);
        zz = xx.^2 + yy.^2;
        
    case 'rosenbrock'
        %-----------------ROSENBROCK FUNCTION 3D PLOT ---------------------
        x = [-50:.1:50];
        y = x;
        k = 25;
        [xx,yy] = meshgrid(x,y);
        zz = ((xx+k)/k-1).^2 + 100*((yy+k)/k-((xx+k)/k).^2).^2;
        
    otherwise
        disp('Invalid objectiveFunction!');
        return
end


% MANIPULATE FITNESSES ----------------------------------------------------
fit = subFitness(swarm,objectiveFunction);
swarmX = swarm(1,:); swarmY = swarm(2,:);

[gBestVal,gBestIndex] = min(fit);
gBest = repmat(swarm(:,gBestIndex),1,noParticle);

pBest = swarm; % initial personal best
pBestFit = subFitness(pBest,objectiveFunction);

BestFitnessEver = zeros(1,maxStep);
MeanFitnessEver = zeros(1,maxStep);




% | UPDATING and TERMINATING SEARCHING PROCESS |||||||||||||||||||||||||||||
for step = 1:maxStep
    % UPDATE VELOCITY /////////////////////////////////////////////////////
    %     % Conventional PSO velocity
    %     velo = (w(step)*velo + c1(step)*rand*(pBest-swarm) + c2(step)*rand*(gBest-swarm));
    % HBPSO velocity
    velo = (w(step)*velo + c1(step)*rand*(pBest-swarm) + c2(step)*rand*(gBest-swarm))...
        + updateBoidVelocity(swarm,velo,staticObs,movingObs);
    
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
    
    
    % UPDATE POSITION OF MOVING OBSTACLES /////////////////////////////////
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
    
    % Find new personal best
    %     pBestFit = subFitness(pBest);
    newFit = fitness(swarm,objectiveFunction); % sensing range is applied
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
    
     
    % CONNECTIVITY CHECK //////////////////////////////////////////////////
    if connectivityCheck == 1
        % Adjacency Matrix
        A = zeros(noParticle,noParticle);
        for i1 = 1:noParticle
            for i2 = (i1+1):noParticle
                if r(swarm(:,i1),swarm(:,i2)) < commuRange
                    A(i1,i2) = 1;
                    A(i2,i1) = 1;
                end
            end
            A(i1,i1) = 0;
        end
        
        % Degree Matrix
        D = zeros(noParticle,noParticle);
        for i1 = 1:noParticle
            for i2 = 1:noParticle
                if (r(swarm(:,i1),swarm(:,i2)) < commuRange) && i1 ~= i2
                    D(i1,i1) = D(i1,i1)+1;
                end
            end
        end
        
        % Laplacian Matrix
        L = D - A;
        
        eigen = sort(eig(L));    lambda2(step) = eigen(2);
    end
     
       
        
        
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
    
    
        
    
    
    % | GRAPHIC |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    % Motion tracking -----------------------------------------------------
    % clear current figure and highlight target
    clf; scatter(0,0,100,red,'filled'); hold on;
    
    % draw obstacles
    scatter(staticObs(1,:),staticObs(2,:),100,black,'filled'); hold on;
    scatter(movingObs(1),movingObs(2),100,black,'filled'); hold on;
    
    xlabel('x coordinator','FontSize',16); ylabel('y coordinator','FontSize',16)
    str = ['Step: ',num2str(step),'\it \color[rgb]{0.2,0.5,0.5} Best Value = ', num2str(gBestVal)];
    title(str,'FontSize',18);
    set(gca,'fontsize',16); hold on;
    
    swarmStep = [swarmX(step,:);swarmY(step,:)]';
    %     viscircles(swarmStep,sensingRange,'EdgeColor','b');
    scatter(swarmX(step,:),swarmY(step,:),90,blue,'filled'); % visual alternation of the above line, draw all the swarm with each particle's size
    
    viscircles(swarmStep,commuRange*sensingRange,'EdgeColor','g','LineWidth',0.1); % draw communication range
    axis([dwLimit_1-1 upLimit_1+1 dwLimit_2-1 upLimit_2+1]); hold on;
    
    if connectivityCheck
        str2 = {['\lambda_2 = ',num2str(lambda2(step))]};
        annotation('textbox', [0.66,0.7,0.1,0.1],'String', str2) % add notation on graph
    end
    
    pause(WAIT/5);
    
    % Get figures for report ----------------------------------------------
    % if step == 1 || step == 5 || step == 20 || step == 40 || step == 60 || step == 100 || step == 150 || step == 200
    %     contour(xx,yy,zz,40);
    %     pause(2);
    % saveas(gcf,num2str(step),'emf'); % save current figure in a defined format
    % end
    
    %     % Trajectory tracking--------------------------------------------
    %     figure;
    %     scatter(0,0,'bo');
    %     xlabel('x coordinator','FontSize',16);
    %     ylabel('ycoordinator','FontSize',16);
    %     a = [num2str(step)];
    %     str = ['Step: ',num2str(step),'\it \color[rgb]{0.2,0.5,0.5} Best Value = ',num2str(gBestVal)];
    %     title(str,'FontSize',12); hold on;
    %     scatter(swarmX(step,:),swarmY(step,:))
    %     axis([dwLimit_1-1 upLimit_1+1 dwLimit_2-1 upLimit_2+1]); hold on;
    %     pause(WAIT/5);   
end
pause(5*WAIT)
hold off

% | DISPLAY RESULTS |||||||||||||||||||||||||||||||||||||||||||||||||||||||
figure;
plot(BestFitnessEver,'color','r','LineWidth',3); hold on; plot(MeanFitnessEver,'LineWidth',3);
grid; title(['\fontsize{18} {\color{blue}Mean of fitness \color{red}vs \color{blue}Best fitness}']);
xlabel('Step ','FontSize',16); ylabel('Fitness','FontSize',16); legend('Best Fitness','MeanFitness');
set(gca,'fontsize',16);
% saveas(gcf,'Searching Performance','emf');
pause(10*WAIT);

if connectivityCheck
    figure;
    plot(lambda2,'LineWidth',2); grid;title('Algebraic Connectivity','Fontsize',20);
    xlabel('Step','FontSize',15); ylabel('\lambda_2','FontSize',15)
    set(gca,'fontsize',16);
    % saveas(gcf,'Algebraic Connectivity','emf');
    % axis([0 maxStep -0.2 noParticle+0.2]);
end

result = ['Best value is ',num2str(gBestVal),', found at (',num2str(gBest(1,1)),',',num2str(gBest(2,1)),')'];
disp(result);
end
% | END ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



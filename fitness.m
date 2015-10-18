function x = fitness(swarm,objectiveFunction)

global noParticle
global range

diagonal = range/sqrt(2);
rg = range*ones(1,noParticle);
dia = diagonal*ones(1,noParticle);
zr = zeros(1,noParticle);

% swarm = dwSpace + rand*spaceSize; % initial position - real position 
swarmN  = swarm + [zr; rg]; % x, y
swarmNW = swarm + [dia; dia];
swarmW  = swarm + [rg; zr];
swarmSW = swarm + [dia; -dia];
swarmS  = swarm + [zr; -rg];
swarmSE = swarm + [-dia; -dia];
swarmE  = swarm + [-rg; zr];
swarmNE = swarm + [-dia; dia];

subN = subFitness(swarmN,objectiveFunction);
subNW = subFitness(swarmNW,objectiveFunction);
subW = subFitness(swarmW,objectiveFunction);
subSW = subFitness(swarmSW,objectiveFunction);
subS = subFitness(swarmS,objectiveFunction);
subSE = subFitness(swarmSE,objectiveFunction);
subE = subFitness(swarmE,objectiveFunction);
subNE = subFitness(swarmNE,objectiveFunction);
sub = subFitness(swarm,objectiveFunction);

tempFitness = [sub; subN; subNW; subW; subSW; subS; subSE; subE; subNE];
x = min(tempFitness);

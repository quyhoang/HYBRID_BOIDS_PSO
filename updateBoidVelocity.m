function veloOutput = updateBoidVelocity(swarm,velo,staticObs,movingObs)

global noParticle
global range

coheAlignSwarm = [swarm,movingObs];
coheAlignSize = size(coheAlignSwarm);

operatingSwarm = [swarm,movingObs,staticObs];
operaSize = size(operatingSwarm);


veloOutput = zeros(operaSize);
mirrorVelo = [velo,zeros(size([staticObs,movingObs]))];
velo = mirrorVelo;
vSeparation = zeros(operaSize);
vAlignment = zeros(operaSize);
vCohesion = zeros(operaSize);
centerCohe = zeros(operaSize); % notice this point
neighbor = 0;

rSeparation = 8*range;
% rSeparation = range + convergence/3;
rCohesion = 20*range;
% rCohesion = 15*range;

% Update vSeparation
for boid1 = 1:operaSize(2)
    for boid2 = 1:operaSize(2)
        if boid2 ~= boid1
            if r(operatingSwarm(:,boid1),operatingSwarm(:,boid2)) < rSeparation
				r12 = r(operatingSwarm(:,boid1),operatingSwarm(:,boid2));
				magnitude = rSeparation - r12;
				connect = operatingSwarm(:,boid1) - operatingSwarm(:,boid2);
				vX = magnitude/r12*connect(1,1);
				vY = magnitude/r12*connect(2,1);
                vSeparation(:,boid1) = vSeparation(:,boid1) + [vX;vY];
            end
        end
    end
end

% Update vCohesion
for boid1 = 1:coheAlignSize(2)
    for boid2 = 1:coheAlignSize(2)
        if boid2 ~= boid1
            if r(operatingSwarm(:,boid1),operatingSwarm(:,boid2)) < rCohesion
                neighbor = neighbor + 1;
                centerCohe(:,boid1) = centerCohe(:,boid1) + operatingSwarm(:,boid2);
            end
        end
    end
    if neighbor ~= 0
        centerCohe(:,boid1) = centerCohe(:,boid1)/neighbor;
        neighbor = 0;
        vCohesion(:,boid1) = (centerCohe(:,boid1) - operatingSwarm(:,boid1))/100; 
        % move 1% towards the center of neighbors
    end
end

neighbor = 0;
% Update vAlignment
for boid1 = 1:coheAlignSize(2)
    for boid2 = 1:coheAlignSize(2)
        if boid2 ~= boid1
            if r(operatingSwarm(:,boid1),operatingSwarm(:,boid2)) < 20*range
                mirrorVelo(:,boid1) = mirrorVelo(:,boid1) + velo(:,boid2);
                neighbor = neighbor + 1; 
            end
        end
    end
    if neighbor ~= 0
        vAlignment(:,boid1) = mirrorVelo(:,boid1)/neighbor;
        neighbor = 0;
        vAlignment(:,boid1) = (vAlignment(:,boid1) - velo(:,boid1))/10;
    end
end

veloOutput = velo + vSeparation*2 + vCohesion*6 + vAlignment*2;
veloOutput = [veloOutput(1,1:noParticle);veloOutput(2,1:noParticle)];
end
% 
% function sqr = sqr(x)
% sqr = x.*x;
% end
% 
% function distance = r(a,b)
% distance = sqrt(sum(sqr(a-b)));
% end


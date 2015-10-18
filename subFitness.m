function sub = subFitness(swarm,objectiveFunction)

global noParticle
sub = ones(1,noParticle);

switch objectiveFunction
    case 'threeHumpCamel'
        %-----------------THREE-HUMP CAMEL FUNCTION -----------------------------
        k = 25;
        for i = 1:noParticle
            sub(i) = 2*(swarm(1,i)/k)^2 - 1.05*(swarm(1,i)/k)^4 + (swarm(1,i)/k)^6/6 + (swarm(1,i)/k)*(swarm(2,i)/k) + (swarm(2,i)/k)^2;
        end
    case 'bohachevsky'
        %-----------------BOHACHEVSKY FUNCTION
        for i = 1:noParticle
            sub(i) = swarm(1,i)^2 + 2*swarm(2,i)^2 - 0.3*cos(3*pi*swarm(1,i) + 4*pi*swarm(2,i)) + 0.3;
        end
    case 'sphere'
        %-----------------SPHERE FUNCTION----------------------------------------
        sub = (sum(swarm.*swarm));
    case 'rosenbrock'
        %-----------------ROSENBROCK FUNCTION -------------------------------------
        k = 25;
        for i = 1:noParticle
            sub(i) = ((swarm(1,i)+k)/k-1).^2 + 100*((swarm(2,i)+k)/k-((swarm(1,i)+k)/k).^2).^2;
        end
    otherwise
        disp('Invalid objectiveFunction!');
        return
end
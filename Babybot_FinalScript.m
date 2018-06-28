clear all; close all; clc

% settings
nbabies = 100;      % number of babies
nbins = 25;     % number of bins
base = 160;      % number of time steps BASELINE
iter = 600;     % number of time steps CONNECTED | default: 450 (15 min * ~30 possible kicks)
exiter = 240;   % number of time steps DISCONNECTED 

% parameters
lr = 0.001;      % learning rate | default: 0.001
d = 0.01;         % decay parameter | default: 0.01
%% START SIMULATING
MCOUNT = zeros(4, nbins); % limbs, bins
for sim = 1:nbabies
    p0 = 1/3*ones(4,3); % initial action probabilities
    x0 = 2*ones(4,1);   % initial state (1 = down, 2 = up, 3 = still)
    p = p0;
    x = x0;
    phist = zeros(4,3,base+iter+exiter);
    xhist = zeros(4,base+iter+exiter);
    ahist = zeros(4,base+iter+exiter);
    %% BASELINE     
    for t = 1:base
        a = zeros(4,1);
        r = rand(4,1);
        a(r<p(:,1)) = 1;
        a(r>=p(:,1)&r<sum(p(:,1:2)')') = 2;
        a(a==0) = 3;
        for limb = 1:4
            p(limb,a(limb)) = p(limb,a(limb))-d*p(limb,a(limb));
            p(limb,:) = p(limb,:)/norm(p(limb,:),1);
        end
       
        % history for the scientist
        ahist(:,t) = a;
        phist(:,:,t) = p;
    end
    %% CONNECTED 
    for t = base+1:base+iter
        % select actions
        a = zeros(4,1);
        r = rand(4,1);
        a(r<p(:,1)) = 1;
        a(r>=p(:,1)&r<sum(p(:,1:2)')') = 2;
        a(a==0) = 3;
        % evaluate effect
        e = zeros(4,1);
        e(([a==1]&[x>=2])|([a==3]&[x<=2])|([a~=2]&[x==2])) = 1;
        % update states
        for limb = find(e')
            if a(limb)==1,x(limb)=x(limb)-1; elseif a(limb)==3,x(limb)=x(limb)+1;end
        end

        % reinforcement OR decay
        if e(1)
            for limb = 1:4
                p(limb,a(limb)) = p(limb,a(limb))+lr*[1-p(limb,a(limb))];
                p(limb,:) = p(limb,:)/norm(p(limb,:),1);
            end
        else
            for limb = 1:4
                p(limb,a(limb)) = p(limb,a(limb))-d*p(limb,a(limb));
                p(limb,:) = p(limb,:)/norm(p(limb,:),1);
            end
        end
        
        % history for the scientist
        ahist(:,t) = a;
        phist(:,:,t) = p;
        xhist(:,t) = x;
    end
    %% DISCONNECTED 
    for t = base+iter+1:base+iter+exiter
        a = zeros(4,1);
        r = rand(4,1);
        a(r<p(:,1)) = 1;
        a(r>=p(:,1)&r<sum(p(:,1:2)')') = 2;
        a(a==0) = 3;
        
        for limb = 1:4
            p(limb,a(limb)) = p(limb,a(limb))-d*p(limb,a(limb));
            p(limb,:) = p(limb,:)/norm(p(limb,:),1);
        end
        
        ahist(:,t) = a;
        phist(:,:,t) = p;
    end
    %% MOVEMENT FREQUENCY
    % bins
    bcountmin = 1;
    bcountmax = 30;
    for bin = 1:nbins
        for limb = 1:4
            mCount(limb, bin) = length(find(ahist(limb, bcountmin:bcountmax) == 1| ahist(limb, bcountmin:bcountmax) == 3));  % per bin count actions (in ahist)
            MCOUNT(limb, bin) = MCOUNT(limb, bin) + mCount(limb, bin);
        end
        bcountmin = bcountmin + 40;
        bcountmax = bcountmax + 40;
    end
end
MCOUNT = MCOUNT/nbabies;
%% PLOT
plot(MCOUNT', 'LineWidth', 5);
xlim([0 nbins])
ylim([18 30])
ylabel('Mean number of kicks')
xlabel('Iteration bins')
legend('left hand (connected)', 'right hand', 'left foot', 'right foot')

save('babybot.mat', 'ahist', 'base', 'd', 'exiter', 'iter', 'lr', 'MCOUNT', 'nbins', 'nbabies', 'phist', 'xhist')

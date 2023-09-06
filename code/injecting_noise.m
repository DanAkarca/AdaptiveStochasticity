%{
Injecting noise
This script runs the generative network model (Vertes et al. 2012, Betzel et
al 2016., Akarca et al. 2021) with heightened stochasticity. The increase
in stochasticity is achieved by setting the probability matrix for novel
connections to zero for a brief period (5% of total connections)
either early, middle, or late in the generative process.

Contact: Dr. Sofia Carozza, sofia.carozza[at]mrc.cbu.cam.ac.uk
%}

function output = injecting_noise(t,...
    parcellation,seed,model,m,el,eu,gl,gu,nsamp,nrep)
%% Define parameters
%{
t: dictates the timing of injecting noise
    1 = early
    2 = middle
    3 = late
parcellation: dictates the costs term
    2 = aal116
seed: dictates the starting conditions
    2 = dHCP rich club 
model: dictates the generative mechanism 
    3 = matching * manuscript
m: dictates the number of connections grown
    400 * manuscript
el, eu: dictates the lower and upper eta limit
    [-4, -0] 
gl, gu: dictates the lower and upper gamma limit
    [0, 1] 
nsamp: dictates the number of evenly spaced samples in the space (must be a ^2 number)
    625 
nrep: dictates the number of repeats undertaken in the space (should be the same as nsamp)
    625 
%}
%% Set directories and paths
% *** set repository directory ***
repo_dir = '/your/path/to/AdaptiveStochasticity';
% change directory to the data directory
cd(strcat(repo_dir,'/data/'));
% load aal116 information
load('aal116_info.mat')
% add brain connectivity toolbox
addpath(genpath('/your/path/to/AdaptiveStochasticity/requirements/'));

%% load neonatal seed
% dHCP seed in AAL
load('dHCP_seed.mat');
%% Load data
% define timing of injecting noise
timing = {'early','middle','late'};
time = timing{t};
switch time
    case 'early' 
        noise_start = round(m*0.05);
        noise_end = round(m*0.1);
    case 'middle'
        noise_start = round(m*0.475);
        noise_end = round(m*0.525);
    case 'late' 
        noise_start = round(m*0.9);
        noise_end = round(m*0.95);
end
%% Initialise model information
% parcellation types
parcellationtypes = string({'dk68','aal116','brainnetome246'});
% model types
modeltypes = string({'sptl',...
    'neighbors','matching',...
    'clu-avg','clu-min','clu-max','clu-diff','clu-prod',...
    'deg-avg','deg-min','deg-max','deg-diff','deg-prod'});
modeltype = modeltypes{model};
% set whether the model is based on powerlaw or exponentials
modelvar = [{'powerlaw'},{'powerlaw'}];
%% Set up parameters
% form an even space based off the eta and gamma limits
x = linspace(el,eu,sqrt(nsamp));
y = linspace(gl,gu,sqrt(nsamp));
params = combvec(x,y)';
% take the cost term
D = squareform(aal116_info.coordinates));
% obtain number of nodes
nnode = size(D,1);
% set seed
if seed == 1
    Aseed = zeros(nnode,nnode);
    mseed = 0;
else if seed == 2
        Aseed = dHCP_seed;
        mseed = nnz(Aseed)/2;
    end
end
% minimum edge value
epsilon     = 1e-5;
% set up output struct
output = struct;
output.parameters = params;
output.networks = zeros(nrep,nsamp,m);
%% Run the grid search: matching
switch modeltype 
    case 'matching'
    % Print text
    disp(sprintf('running the %s model for %g edges within the %s parcellation whilst injecting %s noise...',modeltypes(model),m,parcellationtypes(parcellation),time));
    % Run the model
    for rep = 1:nrep
        for iparam = 1:nsamp
            A = Aseed;
            K = zeros(size(A)); 
            eta = params(iparam,1);
            gam = params(iparam,2);
            K = K + epsilon;
            n = length(D);
            mv1 = modelvar{1};
            mv2 = modelvar{2};
            switch mv1
                case 'powerlaw'
                    Fd = D.^eta;
                case 'exponential'
                    Fd = exp(eta*D);
            end
            switch mv2
                case 'powerlaw'
                    Fk = K.^gam;
                case 'exponential'
                    Fk = exp(gam*K);
            end
            Ff = Fd.*Fk.*~A;
            [u,v] = find(triu(ones(n),1));
            indx = (v - 1)*n + u;
            P = Ff(indx);
            for ii = (mseed+1):(noise_start - 1)
                C = [0; cumsum(P)];
                r = sum(rand*C(end) >= C);
                uu = u(r);
                vv = v(r);

                A(uu,vv) = 1;
                A(vv,uu) = 1;

                updateuu = find(A*A(:,uu));
                updateuu(updateuu == uu) = [];
                updateuu(updateuu == vv) = [];

                updatevv = find(A*A(:,vv));
                updatevv(updatevv == uu) = [];
                updatevv(updatevv == vv) = [];

                c1 = [A(:,uu)', A(uu,:)];
                for i = 1:length(updateuu)
                    j = updateuu(i);
                    c2 = [A(:,j)' A(j,:)];
                    use = ~(~c1&~c2);
                    use(uu) = 0;  use(uu+n) = 0;
                    use(j) = 0;  use(j+n) = 0;
                    ncon = sum(c1(use))+sum(c2(use));
                    if (ncon==0)
                        K(uu,j) = epsilon;
                        K(j,uu) = epsilon;
                    else
                        K(uu,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                        K(j,uu) = K(uu,j);
                    end

                end

                c1 = [A(:,vv)', A(vv,:)];
                for i = 1:length(updatevv)
                    j = updatevv(i);
                    c2 = [A(:,j)' A(j,:)];
                    use = ~(~c1&~c2);
                    use(vv) = 0;  use(vv+n) = 0;
                    use(j) = 0;  use(j+n) = 0;
                    ncon = sum(c1(use))+sum(c2(use));
                    if (ncon==0)
                        K(vv,j) = epsilon;
                        K(j,vv) = epsilon;
                    else
                        K(vv,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                        K(j,vv) = K(vv,j);
                    end
                end
                switch mv2
                    case 'powerlaw'
                        Fk = K.^gam;
                    case 'exponential'
                        Fk = exp(gam*K);
                end
                Ff = Fd.*Fk.*~A;
                P = Ff(indx);
                Ff(isinf(Ff))   = 0; 
            end
            for ii = noise_start:noise_end
                conns = 1:length(P);
                r = randsample(conns(P > 0),1); % choose absolutely randomly
                uu = u(r);
                vv = v(r);

                A(uu,vv) = 1;
                A(vv,uu) = 1;

                updateuu = find(A*A(:,uu));
                updateuu(updateuu == uu) = [];
                updateuu(updateuu == vv) = [];

                updatevv = find(A*A(:,vv));
                updatevv(updatevv == uu) = [];
                updatevv(updatevv == vv) = [];

                c1 = [A(:,uu)', A(uu,:)];
                for i = 1:length(updateuu)
                    j = updateuu(i);
                    c2 = [A(:,j)' A(j,:)];
                    use = ~(~c1&~c2);
                    use(uu) = 0;  use(uu+n) = 0;
                    use(j) = 0;  use(j+n) = 0;
                    ncon = sum(c1(use))+sum(c2(use));
                    if (ncon==0)
                        K(uu,j) = epsilon;
                        K(j,uu) = epsilon;
                    else
                        K(uu,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                        K(j,uu) = K(uu,j);
                    end

                end

                c1 = [A(:,vv)', A(vv,:)];
                for i = 1:length(updatevv)
                    j = updatevv(i);
                    c2 = [A(:,j)' A(j,:)];
                    use = ~(~c1&~c2);
                    use(vv) = 0;  use(vv+n) = 0;
                    use(j) = 0;  use(j+n) = 0;
                    ncon = sum(c1(use))+sum(c2(use));
                    if (ncon==0)
                        K(vv,j) = epsilon;
                        K(j,vv) = epsilon;
                    else
                        K(vv,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                        K(j,vv) = K(vv,j);
                    end
                end
                switch mv2
                    case 'powerlaw'
                        Fk = K.^gam;
                    case 'exponential'
                        Fk = exp(gam*K);
                end
                Ff = Fd.*Fk.*~A;
                P = Ff(indx);
                Ff(isinf(Ff))   = 0;
            end
            for ii = (noise_end+1):m
                C = [0; cumsum(P)];
                r = sum(rand*C(end) >= C);
                uu = u(r);
                vv = v(r);

                A(uu,vv) = 1;
                A(vv,uu) = 1;

                updateuu = find(A*A(:,uu));
                updateuu(updateuu == uu) = [];
                updateuu(updateuu == vv) = [];

                updatevv = find(A*A(:,vv));
                updatevv(updatevv == uu) = [];
                updatevv(updatevv == vv) = [];

                c1 = [A(:,uu)', A(uu,:)];
                for i = 1:length(updateuu)
                    j = updateuu(i);
                    c2 = [A(:,j)' A(j,:)];
                    use = ~(~c1&~c2);
                    use(uu) = 0;  use(uu+n) = 0;
                    use(j) = 0;  use(j+n) = 0;
                    ncon = sum(c1(use))+sum(c2(use));
                    if (ncon==0)
                        K(uu,j) = epsilon;
                        K(j,uu) = epsilon;
                    else
                        K(uu,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                        K(j,uu) = K(uu,j);
                    end

                end

                c1 = [A(:,vv)', A(vv,:)];
                for i = 1:length(updatevv)
                    j = updatevv(i);
                    c2 = [A(:,j)' A(j,:)];
                    use = ~(~c1&~c2);
                    use(vv) = 0;  use(vv+n) = 0;
                    use(j) = 0;  use(j+n) = 0;
                    ncon = sum(c1(use))+sum(c2(use));
                    if (ncon==0)
                        K(vv,j) = epsilon;
                        K(j,vv) = epsilon;
                    else
                        K(vv,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                        K(j,vv) = K(vv,j);
                    end
                end
                switch mv2
                    case 'powerlaw'
                        Fk = K.^gam;
                    case 'exponential'
                        Fk = exp(gam*K);
                end
                Ff = Fd.*Fk.*~A;
                P = Ff(indx);
                Ff(isinf(Ff))   = 0;
            end
            % Keep final network
            b = find(triu(A,1));
            output.networks(rep,iparam,:) = b'; 
        end
    end
    clear eta gam indx A b C r i Fd Ff Fk K u uu v vv iparam x y P ncon 

    case 'sptl'
    % Print text
    disp(sprintf('running the %s model for %g edges within the %s parcellation whilst injecting %s noise...',modeltypes(model),m,parcellationtypes(parcellation),time));
    % Run the model
    for rep = 1:nrep
        for iparam = 1:nsamp
            A = Aseed;
            eta = params(iparam,1);
            n = length(D);
            switch modelvar{1}
                case 'powerlaw'
                    Fd = D.^eta;
                case 'exponential'
                    Fd = exp(eta*D);
            end
            [u,v] = find(triu(ones(n),1));
            indx = (v - 1)*n + u;
            P = Fd(indx).*~A(indx);
            b = zeros(m,1);
            for i = 1:(noise_start-1)
                C = [0; cumsum(P)];
                r = sum(rand*C(end) >= C);
                b(i) = r;
                P = Fd(indx);
                P(b(1:i)) = 0;
            end
            for i = noise_start:noise_end
                conns = 1:length(P);
                r = randsample(conns(P > 0),1); % choose absolutely randomly
                b(i) = r;
                P = Fd(indx);
                P(b(1:i)) = 0;
            end
            for i = (noise_end+1):m
                C = [0; cumsum(P)];
                r = sum(rand*C(end) >= C);
                b(i) = r;
                P = Fd(indx);
                P(b(1:i)) = 0;   
            end
            b = indx(b);
            output.networks(rep,iparam,:) = b'; 
        end
    end
    
    case 'deg-avg'
            % Print text
    disp(sprintf('running the %s model for %g edges within the %s parcellation whilst injecting %s noise...',modeltypes(model),m,parcellationtypes(parcellation),time));
    % Run the model
    for rep = 1:nrep
        for iparam = 1:nsamp
            A = Aseed;
            K = zeros(size(A)); 
            k = sum(A,2);
            eta = params(iparam,1);
            gam = params(iparam,2);
            K = K + epsilon;
            n = length(D);
            [u,v] = find(triu(ones(n),1));
            indx = (v - 1)*n + u;
            Dind = D(indx);
            mv1 = modelvar{1};
            mv2 = modelvar{2};
            switch mv1
                case 'powerlaw'
                    Fd = Dind.^eta;
                case 'exponential'
                    Fd = exp(eta*Dind);
            end
            switch mv2
                case 'powerlaw'
                    Fk = K.^gam;
                case 'exponential'
                    Fk = exp(gam*K);
            end
            P = Fd.*Fk(indx).*~A(indx);
            b = zeros(m,1);
            for i = 1:(noise_start - 1)
                C = [0; cumsum(P)];
                r = sum(rand*C(end) >= C);
                w = [u(r),v(r)];
                k(w) = k(w) + 1;
                switch mv2
                    case 'powerlaw'
                        Fk(:,w) = [((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon].^gam;
                        Fk(w,:) = ([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon].^gam)';
                    case 'exponential'
                        Fk(:,w) = exp([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon]*gam);
                        Fk(w,:) = exp([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon]*gam)';
                end
                P = Fd.*Fk(indx);
                b(i) = r;
                P(b(1:i)) = 0;
            end
            for i = noise_start:noise_end
                conns = 1:length(P);
                r = randsample(conns(P > 0),1); % choose absolutely randomly
                w = [u(r),v(r)];
                k(w) = k(w) + 1;
                switch mv2
                    case 'powerlaw'
                        Fk(:,w) = [((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon].^gam;
                        Fk(w,:) = ([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon].^gam)';
                    case 'exponential'
                        Fk(:,w) = exp([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon]*gam);
                        Fk(w,:) = exp([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon]*gam)';
                end
                P = Fd.*Fk(indx);
                b(i) = r;
                P(b(1:i)) = 0;
            end
            for i = (noise_end+1):m
                C = [0; cumsum(P)];
                r = sum(rand*C(end) >= C);
                w = [u(r),v(r)];
                k(w) = k(w) + 1;
                switch mv2
                    case 'powerlaw'
                        Fk(:,w) = [((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon].^gam;
                        Fk(w,:) = ([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon].^gam)';
                    case 'exponential'
                        Fk(:,w) = exp([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon]*gam);
                        Fk(w,:) = exp([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon]*gam)';
                end
                P = Fd.*Fk(indx);
                b(i) = r;
                P(b(1:i)) = 0;
            end
            % Keep final network
            b = indx(b);
            output.networks(rep,iparam,:) = b'; 
        end
    end
    clear eta gam indx A b C r i Fd Ff Fk K u uu v vv iparam x y P ncon 

    case 'clu-avg'
    % Print text
    disp(sprintf('running the %s model for %g edges within the %s parcellation whilst injecting %s noise...',modeltypes(model),m,parcellationtypes(parcellation),time));
    % Run the model
    for rep = 1:nrep
        for iparam = 1:nsamp
            A = Aseed;
            K = zeros(size(A));
            A = A > 0;
            eta = params(iparam,1);
            gam = params(iparam,2);
            K = K + epsilon;
            n = length(D);
            mv1 = modelvar{1};
            mv2 = modelvar{2};
            switch mv1
                case 'powerlaw'
                    Fd = D.^eta;
                case 'exponential'
                    Fd = exp(eta*D);
            end
            switch mv2
                case 'powerlaw'
                    Fk = K.^gam;
                case 'exponential'
                    Fk = exp(gam*K);
            end
            c = clustering_coef_bu(A);
            k = sum(A,2);
            Ff = Fd.*Fk.*~A;
            [u,v] = find(triu(ones(n),1));
            indx = (v - 1)*n + u;
            P = Ff(indx);
            for ii = 1:(noise_start - 1)
                C = [0; cumsum(P)];
                r = sum(rand*C(end) >= C);
                uu = u(r);
                vv = v(r);
                A(uu,vv) = 1;
                A(vv,uu) = 1;
                k([uu,vv]) = k([uu,vv]) + 1;
                bu = A(uu,:);
                su = A(bu,bu);
                bv = A(vv,:);
                sv = A(bv,bv);
                bth = bu & bv;
                c(bth) = c(bth) + 2./(k(bth).^2 - k(bth));
                c(uu) = nnz(su)/(k(uu)*(k(uu) - 1));
                c(vv) = nnz(sv)/(k(vv)*(k(vv) - 1));
                c(k <= 1) = 0;
                bth([uu,vv]) = true;
                K(:,bth) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')/2 + epsilon;
                K(bth,:) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')'/2 + epsilon;

                switch mv2
                    case 'powerlaw'
                        Ff(bth,:) = Fd(bth,:).*((K(bth,:)).^gam);
                        Ff(:,bth) = Fd(:,bth).*((K(:,bth)).^gam);
                    case 'exponential'
                        Ff(bth,:) = Fd(bth,:).*exp((K(bth,:))*gam);
                        Ff(:,bth) = Fd(:,bth).*exp((K(:,bth))*gam);
                end
                Ff = Ff.*~A;
                P = Ff(indx);
            end
            for ii = noise_start:noise_end
                conns = 1:length(P);
                rng shuffle % set the seed to be shuffled
                r = randsample(conns(P > 0),1); % choose absolutely randomly
                uu = u(r);
                vv = v(r);
                A(uu,vv) = 1;
                A(vv,uu) = 1;
                k([uu,vv]) = k([uu,vv]) + 1;
                bu = A(uu,:);
                su = A(bu,bu);
                bv = A(vv,:);
                sv = A(bv,bv);
                bth = bu & bv;
                c(bth) = c(bth) + 2./(k(bth).^2 - k(bth));
                c(uu) = nnz(su)/(k(uu)*(k(uu) - 1));
                c(vv) = nnz(sv)/(k(vv)*(k(vv) - 1));
                c(k <= 1) = 0;
                bth([uu,vv]) = true;
                K(:,bth) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')/2 + epsilon;
                K(bth,:) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')'/2 + epsilon;

                switch mv2
                    case 'powerlaw'
                        Ff(bth,:) = Fd(bth,:).*((K(bth,:)).^gam);
                        Ff(:,bth) = Fd(:,bth).*((K(:,bth)).^gam);
                    case 'exponential'
                        Ff(bth,:) = Fd(bth,:).*exp((K(bth,:))*gam);
                        Ff(:,bth) = Fd(:,bth).*exp((K(:,bth))*gam);
                end
                Ff = Ff.*~A;
                P = Ff(indx);
            end
            for ii = (noise_end+1):m
                C = [0; cumsum(P)];
                r = sum(rand*C(end) >= C);
                uu = u(r);
                vv = v(r);
                A(uu,vv) = 1;
                A(vv,uu) = 1;
                k([uu,vv]) = k([uu,vv]) + 1;
                bu = A(uu,:);
                su = A(bu,bu);
                bv = A(vv,:);
                sv = A(bv,bv);
                bth = bu & bv;
                c(bth) = c(bth) + 2./(k(bth).^2 - k(bth));
                c(uu) = nnz(su)/(k(uu)*(k(uu) - 1));
                c(vv) = nnz(sv)/(k(vv)*(k(vv) - 1));
                c(k <= 1) = 0;
                bth([uu,vv]) = true;
                K(:,bth) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')/2 + epsilon;
                K(bth,:) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')'/2 + epsilon;

                switch mv2
                    case 'powerlaw'
                        Ff(bth,:) = Fd(bth,:).*((K(bth,:)).^gam);
                        Ff(:,bth) = Fd(:,bth).*((K(:,bth)).^gam);
                    case 'exponential'
                        Ff(bth,:) = Fd(bth,:).*exp((K(bth,:))*gam);
                        Ff(:,bth) = Fd(:,bth).*exp((K(:,bth))*gam);
                end
                Ff = Ff.*~A;
                P = Ff(indx);
           end
            % Keep final network
            b = find(triu(A,1));
            output.networks(rep,iparam,:) = b'; 
        end
    end
    clear eta gam indx A b C r i Fd Ff Fk K u uu v vv x y P ncon 

end

%% Save modelling results 
% Save file
save(sprintf('injecting_noise_%s_%g_%g_%s_%s.mat',parcellationtypes(parcellation),seed,m,modeltypes(model),time),'output','-v7.3');
end
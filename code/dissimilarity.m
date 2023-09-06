%{
Dissimilarity and Robustness Analysis
This script computes the dissimilarity between parameters and repeat runs,
and provides plots. It also runs the robustness analysis.

Contact: Dr. Danyal Akarca, danyal.akarca[at]mrc.cbu.cam.ac.uk
%}
%% add pre-requisites
clear; clc;
% *** set repository directory ***
repo_dir = '/your/path/to/AdaptiveStochasticity';
% change directory to the data directory
cd(strcat(repo_dir,'/data/'));
% load generative models
load('aal116_2_625_625_400_generative_model_matching.mat'); output_aal_seed_400 = output;
% load neonatal seed
load('dHCP_seed.mat');
% load coordinates
load('aal116_info');
coordinates_aal = [aal116.x_mni aal116.y_mni aal116.z_mni]; euclidean_aal = squareform(pdist(coordinates_aal));
% form statistics labels
stat_labels = string({'Clustering','Betweenness','Edge length','Efficiency','Modularity','Communicability'});
% for dissimilarity labels
dissim_labels = string({'\SigmaEuclidean','Consistency'});
% addpath of requirements
addpath(genpath('/Users/da04/Desktop/AdaptiveStochasticity/requirements'));
%% compute node-wise statistics for each network
% Below shows how to compute the node-wise statistics in each network. Alternatively, you can load this pre-computed array
load('global_statistics_six_aal116_seed_400_matching');
%{
% set nnodes
nnode = 116;
% set number of runs
nrun = 625;
% set number of parameters
nparams = 625;
% set number of global statistics measured
ngmeasures = 6;
% set output
output = output_aal_seed_400;
% set euclidean
euclidean = euclidean_aal;
% compute node wise statistics
global_statistics_six_aal116_seed_400_matching = zeros(nrun,nparams,ngmeasures);
% loop over repetitions
for rep = 1:nrun;
    % loop over parameter combinations
    for parametercomb = 1:nparams;
        % form network
        A = zeros(nnode,nnode);
        ind = squeeze(output.networks(rep,parametercomb,:));
        A(ind) = 1;
        A = A + A';
        % compute node statisitcs
        % mean clustering coefficient
        global_statistics_six_aal116_seed_400_matching(rep,parametercomb,1) = mean(clustering_coef_bu(A));
        % mean betweenness centrality
        global_statistics_six_aal116_seed_400_matching(rep,parametercomb,2) = mean(betweenness_bin(A));
        % total edge lengths
        global_statistics_six_aal116_seed_400_matching(rep,parametercomb,3) = sum(euclidean.*A,'all');
        % global efficiency
        global_statistics_six_aal116_seed_400_matching(rep,parametercomb,4) = efficiency_bin(A);
        % modularity q
        [~,global_statistics_six_aal116_seed_400_matching(rep,parametercomb,5)] = modularity_und(A);
        % communicability
        global_statistics_six_aal116_seed_400_matching(rep,parametercomb,6) = mean(expm(A),'all');
        % display
        disp(sprintf('parametercomb %g rep %g node statistics computed',parametercomb,rep));
    end
end
%}
%% compute dissimilarity within parameter combination
% Below shows how to compute within parameter combination dissimilarity. Alternatively, you can load this pre-computed array
load('dissimilarity_aal116_seed_400');
%{
% assign
globalstatistics = global_statistics_six_aal116_seed_400_matching;
output = output_aal_seed_400;
% set nnodes
nnode = 116;
% set number of runs
nrun = 625;
% set number of parameters
nparams = 625;
% set number of connections
nconnect = 400;
% set number of dissimilairty measures computed
ndissim = 2;
% initialise
dissimilarity = zeros(nparams,ndissim);
hdistances = zeros(nparams,nrun,nrun);
cdistances = zeros(nparams,nrun,nrun);
% all combinations
com = combnk(1:nparams,2);
% compute dissimilarity between same parameter runs
for parameterscomb = 1:nparams;
    % get all the runs for this parameter
    data = squeeze(globalstatistics(:,parameterscomb,:)); % check which dimension is parameter, which is repetition
    % normalize the data
    data = normalize(data);
    % calculate the high-dimensional distance
    hdist = squareform(pdist(data));
    % keep
    hdistances(parameterscomb,:,:) = hdist;
    % measure of similarity 1: sum euclidean
    dissimilarity(parameterscomb,1) = sum(hdist,'all');
    % measure of similarity: consistency
    % loop over the runs
    for i = 1:length(com);
        % get the run
        b = squeeze(output.networks(com(i,1),parameterscomb,:));
        c = squeeze(output.networks(com(i,2),parameterscomb,:));
        % compute shared indices
        d = sum(ismember(b,c))/nconnect;
        % compare shared connections
        cdistances(parameterscomb,com(i,1),com(i,2)) = 1-d; % this is now a measure of dissimilarity
    end
    % sum up the shared connectivity across runs
    connections = squeeze(cdistances(parameterscomb,:,:));
    mask = triu(true(size(connections)),1);
    dissimilarity(parameterscomb,2) = mean(connections(mask));
    % display
    disp(sprintf('parametercomb %g dissimilarity computed',parameterscomb));
end
%}
%% produce example networks
% set output
output = output_aal_seed_400;
% set euclidean
euclid = squareform(pdist(coordinates_aal));
% set network
net = 20;
% set run
run = 10;
% set nnode
nnode = 116;
% set coordinates
coordinates = coordinates_aal;
% get network
A = zeros(nnode,nnode);
ind = squeeze(output.networks(net,run,:));
A(ind) = 1;
A = A + A';
g = graph(A);
% plot
% yellow: 246 190 1
% red: 139 0 0 
colplot = [139 0 0]./256;
h = figure; h.Position = [100 100 350 300];
u = plot(g,'XData',coordinates(:,1),...
    'YData',coordinates(:,2),...
    'ZData',coordinates(:,3),...
    'nodelabel',[],...
    'nodecolor',colplot,...
    'markersize',degrees_und(A)+1e-6,...
    'edgecolor',[.7 .7 .7]);
view(270,10); axis off;
% plot the statistical distributions of these simulations
h = figure; h.Position = [100 100 1000 150];
subplot(1,4,1);
histogram(degrees_und(A),'edgecolor','w','facecolor',colplot); 
xlabel('Degree, k'); ylabel('Frequency');
b = gca; b.TickDir = 'out'; b.FontSize = 16; b.FontName = 'Arial'; box off;
subplot(1,4,2);
histogram(clustering_coef_bu(A),'edgecolor','w','facecolor',colplot); 
xlabel('Clustering, c'); ylabel('Frequency');
b = gca; b.TickDir = 'out'; b.FontSize = 16; b.FontName = 'Arial'; box off;
subplot(1,4,3);
histogram(betweenness_bin(A),'edgecolor','w','facecolor',colplot); 
xlabel('Betweenness, b'); ylabel('Frequency');
b = gca; b.TickDir = 'out'; b.FontSize = 16; b.FontName = 'Arial'; box off;
subplot(1,4,4);
histogram(euclid(find(triu(euclid.*A,1))),'edgecolor','w','facecolor',colplot); 
xlabel('Edge length, d'); ylabel('Frequency');
b = gca; b.TickDir = 'out'; b.FontSize = 16; b.FontName = 'Arial'; box off;
%% correlate parameter distance from zero with within-parameter similarity
% set the dissimilarity measure
dissimilarity = dissimilarity_aal116_seed_400;
% *** set dissimiliarity measure to plot *** 
% 1 = Topological dissimilarity, 2 = Embedding dissimiarity
dissim = 1;
% dissimilarity labels
dissim_labels = {'Topological dissimilarity','Embedding dissimilarity'};
% set parameter labels
param_labels = {'\eta','\gamma'};
% get parameter distance from zero
zerodist = [0 0; output.parameters];
zerodist = squareform(pdist(zerodist));
zerodist = zerodist(2:end,1);
% correlate parameter distance from zero with within-parameter dissimilarity
h = figure; h.Position = [100 100 550 450];
x = zerodist; y = dissimilarity(:,dissim);
[r p] = corr(x,y);
scatter(x,y,2000,y,...
    'marker','.');
box off;
% comment in the desired colour scheme
% colormap(brewermap([],"brBG"));
colormap(brewermap([],"BuPu"));
% colormap(brewermap([],"YlGn"));
hold on;
% add a curve
[f gof] = fit(x,y,'poly2');
hold on;
h = plot(f); legend off;
h.Color = [.5 .5 .5]; h.LineWidth = 3;
xlabel('Distance from origin');
ylabel(dissim_labels{dissim});
b = gca; b.TickDir = 'out'; 
b.FontName = 'Arial'; b.FontSize = 16;
% title(sprintf('poly2 R^2=%.3g',100*gof.rsquare));
% place this onto the space
h = figure; h.Position = [100 100 550 450];
scatter(output.parameters(:,1),output.parameters(:,2),2500,dissimilarity(:,dissim),...
    'marker','.'); c = colorbar; c.Label.String = dissim_labels{dissim}; 
% comment in the desired colour scheme
% colormap(brewermap([],"BrBG"));
colormap(brewermap([],"BuPu")); 
% colormap(brewermap([],"YlGn"));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18;
xlabel('\eta'); ylabel('\gamma');
% calculate the variance expalined by each parameter
h = figure; h.Position = [100 100 550 450];
% set parameter
for param = 1:2;
    subplot(1,2,param);
    % visualise
    x = output.parameters(:,param); y = normalize(dissimilarity(:,dissim));
    [r p] = corr(x,y);
    scatter(x,y,2000,y,...
        'marker','.');
    box off;
    % comment in the desired colour scheme
    % colormap(brewermap([],"brBG")); 
    colormap(brewermap([],"BuPu"));
    % colormap(brewermap([],"YlGn"));
    hold on;
    yline(0);
    % add a curve
    [f gof] = fit(x,y,'poly2');
    hold on;
    h = plot(f); legend off;
    h.Color = [.5 .5 .5]; h.LineWidth = 3;
    xlabel(param_labels{param}); % select
    ylabel(dissim_labels{dissim}); 
    b = gca; b.TickDir = 'out'; 
    b.FontName = 'Arial'; b.FontSize = 16;
    % title(sprintf('poly2 R^2=%.3g',100*gof.rsquare));
end
%% plot the specific measures 
% set global statsitics
globalstatistics = global_statistics_six_aal116_seed_400_matching;
% set statistics to plot
x = squeeze(mean(globalstatistics,1));
% number of stats
nstats = 6;
% set statistic to plot
stat = 4;
% set parameter labels
param_labels = {'\eta','\gamma'};
% get parameter distance from zero
zerodist = [0 0; output.parameters];
zerodist = squareform(pdist(zerodist));
zerodist = zerodist(2:end,1);
% correlate parameter distance from zero with the statistics
h = figure; h.Position = [100 100 550 450];
[r p] = corr(zerodist,x(:,stat));
scatter(zerodist,x(:,stat),2000,x(:,stat),...
    'marker','.');
xlabel('Distance from origin');
ylabel(sprintf('%s',stat_labels(stat)));
b = gca; b.TickDir = 'out'; 
b.FontName = 'Arial'; b.FontSize = 14;
box off;
colormap(brewermap([],"BrBG"));
title(sprintf('%r=%.3g, p=%.3g',r,p));
hold on;
% place this onto the space
h = figure; h.Position = [100 100 550 450];
scatter(output.parameters(:,1),output.parameters(:,2),2500,x(:,stat),...
    'marker','.'); c = colorbar; c.Label.String = stat_labels(stat); 
% colormap(brewermap([],"BrBG"));
colormap(brewermap([],"GnBu"));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18;
xlabel('\eta'); ylabel('\gamma');
% plot eta gamma relationships
h = figure; h.Position = [100 100 550 450];
% set parameter
for param = 1:2;
    subplot(1,2,param);
    % visualise
    a = output.parameters(:,param); y = x(:,stat);
    [r p] = corr(a,y);
    scatter(a,y,2000,y,...
        'marker','.');
    box off;
    colormap(brewermap([],"brBG"));
    hold on;
    yline(0);
    % add a curve
    [f gof] = fit(a,y,'poly2');
    hold on;
    h = plot(f); legend off;
    h.Color = [.5 .5 .5]; h.LineWidth = 3;
    xlabel(param_labels{param}); % select
    ylabel(stat_labels(stat)); 
    b = gca; b.TickDir = 'out'; 
    b.FontName = 'Arial'; b.FontSize = 16;
    %title(sprintf('poly2 R^2=%.3g',100*gof.rsquare));
end
% variance explainated by each parameter
concat = [output.parameters x]; 
concat_labels = {'\eta','\gamma','{\itc}','{\itb}','{\itd}','{\ite}','{\itq}','{\itco}'};
% change the order
order = [1 2 3 4 6 5 7 8];
concat = concat(:,order);
concat_labels = concat_labels(order);
% get variance explained
% set model
polymodel = 'poly2';
v = [];
for stat = 1:nstats;
    for param = 1:2;
        [~,gof] = fit(output.parameters(:,param),concat(:,2+stat),polymodel); % assuming order keeps parameters first
        v(stat,param) = gof.rsquare;
    end
end
% get correlation matrix
r = corr(concat);
% visualise
h = figure; h.Position = [100 100 550 450];
imagesc(r); caxis([-1 1]);
colormap(brewermap([],"RdYlBu"));
b = gca; b.TickDir = 'out'; box off; b.FontSize = 14; b.FontName = 'Arial';
xticklabels(concat_labels); 
c = colorbar; c.Label.String = 'r';
% plot eta gamma variance
h = figure; h.Position = [100 100 1000 450];
subplot(1,2,1);
bar(100*v(:,1).^2,.8,'edgecolor','w','facecolor',[.3 .3 .3]);
box off; 
ylabel('\eta R^2 (%)'); xticklabels(concat_labels(3:end));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18;
ylim([0 100]);
subplot(1,2,2);
bar(100*v(:,2).^2,.8,'edgecolor','w','facecolor',[.7 .7 .7]);
box off; 
ylabel('\gamma R^2 (%)'); xticklabels(concat_labels(3:end));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18;
ylim([0 100]);
% plot efficiency versus communicability
% set statistics
eff = x(:,4); comm = x(:,6);
% visualise
h = figure; h.Position = [100 100 550 450];
[r p] = corr(eff,comm);
scatter(eff,comm,1500,...
    'marker','.',...
    'markeredgecolor',[.5 .5 .5]);
xlabel('Efficiency');
ylabel('Communicability');
b = gca; b.TickDir = 'out'; 
b.FontName = 'Arial'; b.FontSize = 16;
box off;
colormap(brewermap([],"BrBG"));
title(sprintf('%r=%.3g, p=%.3g',r,p));
%% compute network robustness across the parameter space
% Below shows the robustness analysis. Alternatively, you can load this pre-computed array
load('robustness_aal116_seed_400_matching');
%{
% set nnodes
nnode = 116;
% set number of runs
nrun = 625;
% set number of parameters
nparams = 625;
% set output
output = output_aal_seed_400;
% set euclidean
euclidean = euclidean_aal;
% set the number of nodes to remove
nremove = 0.25*nnode;
% set number of robustness regimes to test
nregime = 3;
% compute node wise statistics
robustness = zeros(nrun,nparams,nremove,nregime);
robustness_summary = zeros(nrun,nparams,nregime);
% robustness operates by removing high degree nodes according to a regime
% and then calculating the change of communication or efficiency
% loop over repetitions
for rep = 1:nrun;
    % loop over parameter combinations
    for parametercomb = 1:nparams;
        % loop over regimes
        for regime = 1:nregime;
            % initialise
            b = [];
             % form network
            A = zeros(nnode,nnode);
            ind = output.networks(rep,parametercomb,:);
            A(ind) = 1;
            A = A + A';
            % targeted attack (high to low)
            if regime == 1;
                % calculate the degree distribution
                k = degrees_und(A);
                % rank them
                [~,i] = sort(k,'descend');
                % remove in order and recalculate
                for node = 1:nremove;
                    % remove in order of most communicable to least
                    A(i(node),:) = 0;
                    A(:,i(node)) = 0;
                    b(node) = sum(expm(A),'all');
                end
            end
            % targeted attack (low to high)
            if regime == 2;
                % calculate the degree distribution
                k = degrees_und(A);
                % rank them
                [~,i] = sort(k,'ascend');
                % remove in order and recalculate
                for node = 1:nremove;
                    % remove in order of most communicable to least
                    A(i(node),:) = 0;
                    A(:,i(node)) = 0;
                    b(node) = sum(expm(A),'all');
                end
            end
            % random attack
            if regime == 3;
                % get a random order of nodes
                i = randperm(nnode);
                % remove in order and recalculate
                for node = 1:nremove;
                    % remove in order of most communicable to least
                    A(i(node),:) = 0;
                    A(:,i(node)) = 0;
                    b(node) = sum(expm(A),'all');
                end
            end
            % fit a model and keep the gradient
            % keep the data
            robustness(rep,parametercomb,:,regime) = b;
            model = fit([1:nremove]',log(b'),'poly1'); % natural logarithm
            robustness_bet(rep,parametercomb,regime) = model.p1; % gradient
            robustness_abs(rep,parametercomb,regime) = b(1)-b(end); % absolute drop
            robustness_rel(rep,parametercomb,regime) = 1-(b(end)/b(1)); % relative drop
            robustness_end(rep,parametercomb,regime) = b(end); % retained
        end
        % display
        disp(sprintf('parametercomb %g rep %g robustness regimes computed',parametercomb,rep));
    end
end
%}
%% examine robustness
% set data
robustness_data = robustness_aal116_seed_400_matching;
% put data together
robustness = robustness_data.robustness;
beta = robustness_data.beta;
absolute = robustness_data.absolute;
relative = robustness_data.relative;
ending = robustness_data.end;
robust_dataplot = {beta,absolute,relative,ending};
% set parameter labels
param_labels = {'\eta','\gamma'};
% get parameter distance from zero
zerodist = [0 0; output.parameters];
zerodist = squareform(pdist(zerodist));
zerodist = zerodist(2:end,1);
% set regime labels
regime_labels = string({'Targetted (H2L) attacks','Targetted (L2H) attacks','Random attacks'});
% *** set regime ***
regime = 3;
% set data labels
data_labels = string({'\beta coefficient','{\delta\itC_{i,j}} loss','% {\itC_{i,j}} loss','Retained {\itC_{i,j}}'});
% *** set robustness data to plot ***
dataplot = 1; 
% get data and take the derivative measure
robust_regime = squeeze(mean(robust_dataplot{dataplot}(:,:,regime),1))';
% correlate parameter distance from zero with the robustness
h = figure; h.Position = [100 100 550 450];
[r p] = corr(zerodist,robust_regime);
scatter(zerodist,robust_regime,2000,robust_regime,...
    'marker','.');
% sgtitle(sprintf('%s',regime_labels(regime)));
xlabel('Distance from origin'); ylabel(data_labels(dataplot));
b = gca; b.TickDir = 'out'; 
b.FontName = 'Arial'; b.FontSize = 14;
box off;
colormap(brewermap([],"RdYlBu"));
title(sprintf('%r=%.3g, p=%.3g',r,p));
hold on;
% place this onto the space
h = figure; h.Position = [100 100 550 450];
scatter(output.parameters(:,1),output.parameters(:,2),2500,robust_regime,...
    'marker','.'); c = colorbar; c.Label.String = data_labels(dataplot); colormap(brewermap([],"RdYlBu"));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18; 
% set a manual caxis if wanted
% caxis([-0.65 -0.1]);
% sgtitle(regime_labels(regime));
xlabel('\eta'); ylabel('\gamma');
% plot eta gamma relationships
h = figure; h.Position = [100 100 550 450];
% set parameter
for param = 1:2;
    subplot(1,2,param);
    % visualise
    a = output.parameters(:,param); y = robust_regime;
    [r p] = corr(a,y);
    scatter(a,y,2000,y,...
        'marker','.');
    box off;
    colormap(brewermap([],"RdYlBu"));
    hold on;
    % add a curve
    [f gof] = fit(a,y,'poly2');
    hold on;
    h = plot(f); legend off;
    h.Color = [.5 .5 .5]; h.LineWidth = 3;
    xlabel(param_labels{param});
    ylabel(data_labels(dataplot));
    b = gca; b.TickDir = 'out'; 
    b.FontName = 'Arial'; b.FontSize = 16;
    % title(sprintf('poly2 R^2=%.3g',100*gof.rsquare));
    % sgtitle(sprintf('%s',regime_labels(regime)));
end
% plot the most and least robust trajectories
[h hi] = max(robust_regime); x = squeeze(robustness(:,hi,:,regime)); x = log(x'); % natural log
[l li] = min(robust_regime); y = squeeze(robustness(:,li,:,regime)); y = log(y'); % natural log
% get limits for the plot
u = [max(x,[],'all'),max(y,[],'all')];
v = [min(x,[],'all'),min(y,[],'all')];
climits = [min(v) max(u)];
% place data together for the loop
robust_data = {x y};
% get colours
z = brewermap([],"RdYlBu");
col(1,:) = z(end,:);
col(2,:) = z(1,:);
% get indexes
in = [hi li];
% plot high and low robustness
% plot low robustness first then high
plotorder = [2 1];
h = figure; h.Position = [100 100 550 450];
for ix = 1:2;
    subplot(1,2,ix);
    % visualise
    stdshade(robust_data{plotorder(ix)}',.6,col(plotorder(ix),:));
    box off;
    colormap(brewermap([],"RdYlBu"));
    % title(sprintf('/eta=%g, /gamma=%g',output.parameters(in(param),1),output.parameters(in(param),2)));
    xlabel('Removed node');
    ylabel('ln({\itC_{i,j}})');
    % sgtitle(sprintf('%s',regime_labels(regime)));
    b = gca; b.TickDir = 'out'; 
    b.FontName = 'Arial'; b.FontSize = 16;
    ylim(climits);
end
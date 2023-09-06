%{
Noise Analysis
This script analyses the output of injecting_noise, which artificially
heightened the stochasticity of the generative model either early, middle, 
or late in the generative process. It assesses the impact of injected
noise on the resulting topological dissimilarity and embedding
dissimilarity across multiple runs of the same parameter combinations,
and creates figures to visualise these results.

Contact: Dr. Sofia Carozza, sofia.carozza[at]mrc.cbu.cam.ac.uk
%}
%% set parameters
parcellation = 'aal116'; 
m = 400;
model = 'matching';
nrun = 625;
nparams = 625;
seed = 2; % 2 for neonatal rich club seed
%% set paths, add tools and import data
% *** set repository directory ***
repo_dir = '/your/path/to/AdaptiveStochasticity';
% change directory to the data directory
cd(strcat(repo_dir,'/data/'));
% add bct and color brewer to path
addpath(genpath('/your/path/to/AdaptiveStochasticity/requirements/'));
% load results from injecting_noise script
outputs = cell(1,3);
load(strcat('injecting_noise_',parcellation,'_',num2str(seed),'_',num2str(m),'_',model,'_early')); 
outputs{1} = output; clear output
load(strcat('injecting_noise_',parcellation,'_',num2str(seed),'_',num2str(m),'_',model,'_middle')); 
outputs{2} = output; clear output
load(strcat('injecting_noise_',parcellation,'_',num2str(seed),'_',num2str(m),'_',model,'_late')); 
outputs{3} = output; 
% load distance matrix
load('red_qsiprep.mat');
euclidean = squareform(pdist(red_qsiprep.aal116.coordinates));
clear red_qsiprep

%% compute global statistics for each network
% set stats
nnode = size(euclidean,1);
% set number of global statistics measured
ngmeasures = 5;
% create cell array to save statistics
global_statistics = cell(1,3);
% loop over noise stages
for t = 1:3
    statistics = zeros(625,625,ngmeasures);
    output = outputs{t};
    % loop over repetitions
    for rep = 1:nrun
        % loop over parameter combinations
        for parametercomb = 1:nparams
            % form network
            A = zeros(nnode,nnode);
            ind = output.networks(rep,parametercomb,:);
            A(ind) = 1;
            A = A + A';
            % compute node statisitcs
            % mean clustering coefficient
            statistics(rep,parametercomb,1) = mean(clustering_coef_bu(A));
            % mean betweenness centrality
            statistics(rep,parametercomb,2) = mean(betweenness_bin(A));
            % total edge lengths
            statistics(rep,parametercomb,3) = sum(sum(euclidean.*A));
            % global efficiency
            statistics(rep,parametercomb,4) = efficiency_bin(A);
            % modularity q
            [~,statistics(rep,parametercomb,5)] = modularity_und(A);
            % display
            disp(sprintf('parametercomb %g rep %g node statistics computed',parametercomb,rep));
        end
    end
    global_statistics{t} = statistics;
end
save(strcat('global_statistics_',parcellation,'_',num2str(seed),'_',num2str(m),'_',model,'.mat'),'global_statistics');
clear A output t rep nnodes statistics ngmeasures parametercomb
%% compute dissimilarity within parameter combination
% load data if previously run
% load(strcat(datapath,'global_statistics_',parcellation,'_',num2str(m),'_',model,'.mat'));
% otherwise form statistics labels
stat_labels = string({'Clustering','Betweenness','Edge length','Efficiency','Modularity'});
% set number of dissimilarity measures computed
ndissim = 2;
% write the label
dissim_labels = string({'\SigmaEuclidean','Consistency'});
% loop through all three noise points
dissimilarities = cell(1,3);
for t = 1:3
    % get data for that time point
    statistics = global_statistics{t};
    output = outputs{t};
    % initialise variables
    dissimilarity = zeros(nparams,ndissim);
    hdistances = zeros(nparams,nrun,nrun);
    cdistances = zeros(nparams,nrun,nrun);
    % all combinations
    com = combnk(1:nrun,2);
    % compute dissimilarity between same parameter runs
    for parameterscomb = 1:nparams
        % get all the runs for this parameter
        data = squeeze(statistics(:,parameterscomb,:));
        % normalize the data
        data = normalize(data);
        % calculate the high-dimensional distance
        hdist = squareform(pdist(data));
        % keep high-dimensional distance
        hdistances(parameterscomb,:,:) = hdist;
        % measure of dissimilarity 1: sum euclidean
        dissimilarity(parameterscomb,1) = sum(sum(hdist));
        % measure of dissimilarity 2: consistency
        % loop over the runs
        for i = 1:length(com)
            % get the run
            b = squeeze(output.networks(com(i,1),parameterscomb,:));
            c = squeeze(output.networks(com(i,2),parameterscomb,:));
            % compute shared indices
            d = sum(ismember(b,c))/m;
            % compare shared connections
            cdistances(parameterscomb,com(i,1),com(i,2)) = d; % note this is actually a measure of similarity
        end
        % sum up the shared connectivity across runs
        connections = squeeze(cdistances(parameterscomb,:,:));
        mask = triu(true(size(connections)),1);
        dissimilarity(parameterscomb,2) = 1 - mean(connections(mask));
        % display
        disp(sprintf('parametercomb %g dissimilarity computed for time %g',parameterscomb,t));
    end
    dissimilarities{t} = dissimilarity;
end
save(strcat('dissim__',parcellation,'_',num2str(seed),'_',num2str(m),'_',model,'.mat'),'dissimilarities');
clear t parameterscomb data hdist ndissim connections mask hdistancesi b c d cdistances dissimilarity 
%% correlate parameter distance from zero with within-parameter similarity
% load control 
load(strcat('dissimilarity_',parcellation,'_seed_',num2str(m)));
dissimilarity_control = dissimilarity_aal116_seed_400; clear dissimilarity_aal116_seed_400

% extract dissimilarities for each time point 
dissimilarity_early = dissimilarities{1}; dissimilarity_middle = dissimilarities{2}; dissimilarity_late = dissimilarities{3};
%% set dissimiliarity measure to plot (1 for embedding, 2 for topology)
dissim = 2; 
if dissim == 2
    colormap = "BuPu";
else 
    colormap = "YlGn";
end

%% generate dissimilarity scatterplots
output_params = [2 1];
param_labels = {'\gamma','\eta'};
figurepath = strcat(repo_dir,'/',parcellation,'_',num2str(seed),'_',num2str(m),'_',model,'_',num2str(dissim),'/');
% get parameter distance from zero
zerodist = [0 0; output.parameters];
zerodist = squareform(pdist(zerodist));
zerodist = zerodist(2:end,1);
% correlate parameter distance from zero with within-parameter dissimilarity
h = figure; 
h.Position = [10 100 350 400];
hold on
x = output.parameters(:,output_params(dissim)); y = dissimilarity_control(:,dissim);
scatter(x,y,400,y,...
    'marker','.');
box off;
ylim([0.99*min(dissimilarity_control(:,dissim)) 1.01*max(dissimilarity_control(:,dissim))])
colormap(brewermap([],colormap));
caxis([min(dissimilarity_control(:,dissim)) max(dissimilarity_control(:,dissim))]);
xlim([min(output.parameters(:,output_params(dissim))) max(output.parameters(:,output_params(dissim)))])
[f gof] = fit(x,y,'poly2');
title({'No noise',sprintf('R^2=%.3g',100*gof.rsquare)});
h = plot(f); legend off;
h.Color = [.5 .5 .5]; h.LineWidth = 3;
xlabel(param_labels{dissim});
ylabel('Dissimilarity');
b = gca; b.TickDir = 'out'; 
b.FontName = 'Arial'; b.FontSize = 16;
hold off
print(strcat(figurepath,'dissimilarity_scatter_control',num2str(dissim)),'-djpeg','-r300');

h = figure; 
h.Position = [10 100 1200 400];
subplot(1,3,1)
hold on
x = output.parameters(:,output_params(dissim)); y = dissimilarity_early(:,dissim);
scatter(x,y,400,y,'marker','.');
box off;
colormap(brewermap([],colormap));
caxis([min(dissimilarity_control(:,dissim)) max(dissimilarity_control(:,dissim))]);
xlim([min(output.parameters(:,output_params(dissim))) max(output.parameters(:,output_params(dissim)))])
ylim([0.99*min(dissimilarity_control(:,dissim)) 1.01*max(dissimilarity_control(:,dissim))])
% add a curve
[f gof] = fit(x,y,'poly2');
title({'Early noise',sprintf('R^2=%.3g',100*gof.rsquare)});
h = plot(f); legend off;
h.Color = [.5 .5 .5]; h.LineWidth = 3;
xlabel(param_labels{dissim});
ylabel('Dissimilarity');
b = gca; b.TickDir = 'out'; 
b.FontName = 'Arial'; b.FontSize = 16;
hold off

subplot(1,3,2)
hold on
x = output.parameters(:,output_params(dissim)); y = dissimilarity_middle(:,dissim);
scatter(x,y,400,y,...
    'marker','.');
box off;
ylim([0.99*min(dissimilarity_control(:,dissim)) 1.01*max(dissimilarity_control(:,dissim))])
colormap(brewermap([],colormap));
caxis([min(dissimilarity_control(:,dissim)) max(dissimilarity_control(:,dissim))]);
xlim([min(output.parameters(:,output_params(dissim))) max(output.parameters(:,output_params(dissim)))])
[f gof] = fit(x,y,'poly2');
h = plot(f); legend off;
title({'Middle noise',sprintf('R^2=%.3g',100*gof.rsquare)});
h.Color = [.5 .5 .5]; h.LineWidth = 3;
xlabel(param_labels{dissim});
ylabel('Dissimilarity');
b = gca; b.TickDir = 'out'; 
b.FontName = 'Arial'; b.FontSize = 16;
hold off

subplot(1,3,3)
hold on
x = output.parameters(:,output_params(dissim)); y = dissimilarity_late(:,dissim);
scatter(x,y,400,y,...
    'marker','.');
box off;
ylim([0.99*min(dissimilarity_control(:,dissim)) 1.01*max(dissimilarity_control(:,dissim))])
colormap(brewermap([],colormap));
caxis([min(dissimilarity_control(:,dissim)) max(dissimilarity_control(:,dissim))]);
xlim([min(output.parameters(:,output_params(dissim))) max(output.parameters(:,output_params(dissim)))])
[f gof] = fit(x,y,'poly2');
title({'Late noise',sprintf('R^2=%.3g',100*gof.rsquare)});
h = plot(f); legend off;
h.Color = [.5 .5 .5]; h.LineWidth = 3;
xlabel(param_labels{dissim});
ylabel('Dissimilarity');
b = gca; b.TickDir = 'out'; 
b.FontName = 'Arial'; b.FontSize = 16;
hold off
print(strcat(figurepath,'/dissimilarity_scatter',num2str(dissim)),'-djpeg','-r300');
%% generate heatmaps of parameter space and dissimilarity
h = figure;
h.Position = [10 100 400 340];
hold on
scatter(output.parameters(:,1),output.parameters(:,2),1000,dissimilarity_control(:,dissim),...
    'marker','.'); colormap(brewermap([],colormap)); c = colorbar; %c.Label.String = 'Dissimilarity'; 
caxis([min(dissimilarity_control(:,dissim)) max(dissimilarity_control(:,dissim))]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18;
xlabel('\eta'); ylabel('\gamma'); xticklabels([-4 -2 0]); 
yticks([0 1]); yticklabels([0 1]);
hold off
print(strcat(figurepath,'dissimilarity_heatmap_control'),'-djpeg','-r300');

h = figure;
h.Position = [10 100 1200 340];
subplot(1,3,1)
hold on
scatter(output.parameters(:,1),output.parameters(:,2),1000,dissimilarity_early(:,dissim),...
    'marker','.'); colormap(brewermap([],colormap)); 
caxis([min(dissimilarity_control(:,dissim)) max(dissimilarity_control(:,dissim))]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18;
xlabel('\eta'); ylabel('\gamma'); xticklabels([-4 -2 0]); 
yticks([0 1]); yticklabels([0 1]);
hold off
subplot(1,3,2)
hold on
scatter(output.parameters(:,1),output.parameters(:,2),1000,dissimilarity_middle(:,dissim),...
    'marker','.'); colormap(brewermap([],colormap)); 
caxis([min(dissimilarity_control(:,dissim)) max(dissimilarity_control(:,dissim))]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18;
xlabel('\eta'); ylabel('\gamma'); xticklabels([-4 -2 0]); 
yticks([0 1]); yticklabels([0 1]);
hold off
subplot(1,3,3)
hold on
scatter(output.parameters(:,1),output.parameters(:,2),1000,dissimilarity_late(:,dissim),...
    'marker','.'); colormap(brewermap([],colormap)); 
caxis([min(dissimilarity_control(:,dissim)) max(dissimilarity_control(:,dissim))]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18;
xlabel('\eta'); ylabel('\gamma'); xticklabels([-4 -2 0]); 
yticks([0 1]); yticklabels([0 1]);
hold off

print(strcat(figurepath,'dissimilarity_heatmap'),'-djpeg','-r300');
%% heatmap of dissimilarity relative to control
h = figure;
h.Position = [10 100 800 340];
hold on
scatter(output.parameters(:,1),output.parameters(:,2),1000,...
    (dissimilarity_early(:,dissim)-dissimilarity_control(:,dissim)),...
    'marker','.'); colormap(brewermap([],colormap)); 
c = colorbar; c.Label.String = '\Delta Dissimilarity'; 
caxis([min(dissimilarity_early(:,dissim)-dissimilarity_control(:,dissim)) max(dissimilarity_early(:,dissim)-dissimilarity_control(:,dissim))]);...
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18;
xlabel('\eta'); ylabel('\gamma'); xticklabels([-4 -2 0]); 
yticks([0 1]); yticklabels([0 1]);
hold off
print(strcat(figurepath,'dissimilarity_heatmap_relative_colorbar'),'-djpeg','-r300');

h = figure;
h.Position = [10 100 1200 340];
subplot(1,3,1)
hold on
scatter(output.parameters(:,1),output.parameters(:,2),1000,...
    (dissimilarity_early(:,dissim)-dissimilarity_control(:,dissim)),...
    'marker','.'); colormap(brewermap([],colormap)); 
caxis([min(dissimilarity_early(:,dissim)-dissimilarity_control(:,dissim)) max(dissimilarity_early(:,dissim)-dissimilarity_control(:,dissim))]);...
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18;
xlabel('\eta'); ylabel('\gamma'); xticklabels([-4 -2 0]); 
yticks([0 1]); yticklabels([0 1]);
hold off
subplot(1,3,2)
hold on
scatter(output.parameters(:,1),output.parameters(:,2),1000,...
    (dissimilarity_middle(:,dissim)-dissimilarity_control(:,dissim)),...
    'marker','.'); colormap(brewermap([],colormap)); 
caxis([min(dissimilarity_early(:,dissim)-dissimilarity_control(:,dissim)) max(dissimilarity_early(:,dissim)-dissimilarity_control(:,dissim))]);...
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18;
xlabel('\eta'); ylabel('\gamma'); xticklabels([-4 -2 0]); 
yticks([0 1]); yticklabels([0 1]);
hold off
subplot(1,3,3)
hold on
scatter(output.parameters(:,1),output.parameters(:,2),1000,...
    (dissimilarity_late(:,dissim)-dissimilarity_control(:,dissim)),...
    'marker','.'); colormap(brewermap([],colormap)); 
caxis([min(dissimilarity_early(:,dissim)-dissimilarity_control(:,dissim)) max(dissimilarity_early(:,dissim)-dissimilarity_control(:,dissim))]);...
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 18;
xlabel('\eta'); ylabel('\gamma'); xticklabels([-4 -2 0]); 
yticks([0 1]); yticklabels([0 1]);
hold off

print(strcat(figurepath,'dissimilarity_heatmap_relative'),'-djpeg','-r300');
%% scatter plots of delta dissimilarity and parameters
figurepath = strcat(repo_dir,'/',parcellation,'_',num2str(seed),'_',num2str(m),'_',model,'_',num2str(dissim),'/');
delta_dissim = [dissimilarity_early(:,dissim)-dissimilarity_control(:,dissim),...
    dissimilarity_middle(:,dissim)-dissimilarity_control(:,dissim),...
    dissimilarity_late(:,dissim)-dissimilarity_control(:,dissim)];
[p,t,stats] = anova1(delta_dissim);

% relative dissimilarity
h = figure;
h.Position = [10 100 900 250];
subplot(1,3,1)
if dissim == 2
    set(gca,'xdir','reverse')
else
    set(gca,'xdir','normal')
end
hold on
x = output.parameters(:,output_params(dissim));
y = delta_dissim(:,1);
scatter(x,y,30,y)
[f gof] = fit(x,y,'poly2');
h = plot(f); legend off;
h.Color = [.5 .5 .5]; h.LineWidth = 3;
xlabel(param_labels{dissim});
title({'Early noise',sprintf('R^2=%.3g%%',100*gof.rsquare)});
ylabel('\Delta Dissimilarity'); 
colormap(brewermap([],"YlGn"));
caxis([0 max(delta_dissim(:))])
ylim([0 max(delta_dissim(:))])
hold off
subplot(1,3,2)
if dissim == 2
    set(gca,'xdir','reverse')
else
    set(gca,'xdir','normal')
end
hold on
x = output.parameters(:,output_params(dissim));
y = delta_dissim(:,2);
scatter(x,y,30,y)
[f gof] = fit(x,y,'poly2');
h = plot(f); legend off;
h.Color = [.5 .5 .5]; h.LineWidth = 3;
xlabel(param_labels{dissim});
title({'Middle noise',sprintf('R^2=%.3g%%',100*gof.rsquare)});
ylabel('\Delta Dissimilarity'); 
colormap(brewermap([],"YlGn"));
caxis([0 max(delta_dissim(:))])
ylim([0 max(delta_dissim(:))])
hold off
subplot(1,3,3)
if dissim == 2
    set(gca,'xdir','reverse')
else
    set(gca,'xdir','normal')
end
hold on
x = output.parameters(:,output_params(dissim));
y = delta_dissim(:,3);
scatter(x,y,30,y)
[f gof] = fit(x,y,'poly2');
h = plot(f); legend off;
h.Color = [.5 .5 .5]; h.LineWidth = 3;
xlabel(param_labels{dissim});
title({'Late noise',sprintf('R^2=%.3g%%',100*gof.rsquare)});
ylabel('\Delta Dissimilarity'); 
colormap(brewermap([],"YlGn"));
caxis([0 max(delta_dissim(:))])
ylim([0 max(delta_dissim(:))])
hold off
print(strcat(figurepath,'delta_dissim'),'-djpeg','-r300');



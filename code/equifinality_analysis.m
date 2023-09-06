%{
Equifinality Analysis
This script analyses the output of run_missclass. It tests two hypothesized
contributors to equifinality: distance between the gamma parameters of the
two sets of simulations, and the inherent stochasticity of the model. It 
produces plots to convey these results.

Contact: Dr. Sofia Carozza, sofia.carozza[at]mrc.cbu.cam.ac.uk
%}

%% set parameters and paths
% set parameters
parcellation = 'aal116';
m = 400;
model = 'matching';
nrun = 625;
nparams = 625;
seed = 2;
dissim = 1;
% *** set repository directory ***
repo_dir = '/your/path/to/AdaptiveStochasticity';
% change directory to the data directory
cd(strcat(repo_dir,'/data/'));
% add bct and color brewer to path
addpath(genpath('/your/path/to/AdaptiveStochasticity/requirements/'));
% set figure path
figurepath = strcat(repo_dir,'/',parcellation,'_',num2str(seed),'_',num2str(m),'_',model,'_',num2str(dissim),'/');
%% load data
% load generative models
load('aal116_2_625_625_400_generative_model_matching.mat'); 
load('aal116_info');
coordinates = [aal116.x_mni aal116.y_mni aal116.z_mni]; euclidean = squareform(pdist(coordinates));
load('dissimilarity_aal116_seed_400.mat');
load('global_statistics_aal116_seed_400_matching.mat');
global_statistics = global_statistics_aal116_seed_400_matching; clear global_statistics_aal116_seed_400_matching
dissimilarity = dissimilarity_aal116_seed_400; clear dissimilarity_aal116_seed_400 aal116

% form statistics labels
stat_labels = string({'Clustering','Betweenness','Edge length','Efficiency','Modularity'});
% for dissimilarity labels
dissim_labels = string({'\SigmaEuclidean','Consistency'});
%% load missclassification values for all possible pairs
missclass_625 = zeros(nparams,nparams);
for p = 1:nparams
    load(strcat('missclass_625/missclass_',num2str(p),'.mat'));
    missclass_625(p,:) = missclass;
    clear missclass
end
clear p
%% compute parameter distance and total topological dissimilarity
param_dist = zeros(nparams,nparams);
gamma_dist = zeros(nparams,nparams);

load('625_params.mat');
for p = 1:nparams
    for n = 1:nparams
        param_dist(p,n) = pdist([parameters(p,:);...
            parameters(n,:)],'euclidean');
        gamma_dist(p,n) = abs(parameters(p,2) - parameters(n,2));
    end
end
clear p n

%% plot success as a function of dissimilarity
missclass_625_0 = missclass_625;
missclass_625_0(isnan(missclass_625_0)) = 0;

figure; hold on
x = normalize(dissimilarity(:,1),'range',[0,100]);
y = mean(missclass_625_0,2);
scatter(x,y,'MarkerFaceColor',[0 .4 .7])
[f, gof] = fit(x,y,'poly2');
disp(f)
h = plot(f); legend off;
h.Color = [.5 .5 .5]; h.LineWidth = 3;
ylabel('Mean misclassification rate','FontSize',14)
xlabel('Topological dissimilarity (% of max)','FontSize',14)
title({'Equifinality by simulation stochasticity',sprintf('R^2=%.3g%%',100*gof.rsquare)},...
    'FontSize',16);
hold off
print(strcat(figurepath,'missclass_dissimilarity'),'-djpeg','-r300');

%% equifinality by gamma distance
hold on
figure
x = reshape(gamma_dist,1,[])';
y = reshape(missclass_625,1,[])';
densityplot(x,y,[5,5])
title({'Equifinality by \gamma distance'},...
    'FontSize',16);
xlabel('Difference between \gamma values','FontSize',14)
ylabel('Misclassification rate','FontSize',14)
colormap parula
hold off
print(strcat(figurepath,'missclass_distance'),'-djpeg','-r300');

%% plot heatmap of classifier success
figure;
scatter(parameters(:,1),parameters(:,2),1000,mean(missclass_625_0,2),...
    'marker','.'); colormap(brewermap([],"PuBu")); c = colorbar; c.Label.String = 'Mean misclassification rate'; 
caxis([min(mean(missclass_625_0,2)) max(mean(missclass_625_0,2))]);
c.Location = 'westoutside';c.FontSize =14;
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('\eta'); ylabel('\gamma');xticks([-4 -2 0]); 
yticks([0 1]); yticklabels([0 1]);
print(strcat(figurepath,'missclass_heatmap'),'-djpeg','-r300');

%% compute r2 values
[r p] = corrcoef(normalize(dissimilarity(:,1),'range',[0,100]),mean(missclass_625_0,2))
[r p] = corrcoef(reshape(gamma_dist,1,[])',reshape(missclass_625_0,1,[])')

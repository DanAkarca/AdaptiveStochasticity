%{
Run MissClass
This script runs a Support Vector Machine (SVM) on the topological
characteristics of sets of repeated runs of the
generative network model. Each set was obtained with different wiring 
parameters (eta and gamma). As such, the success of the SVM indicates the 
similarity, across sets, in global topology of the networks. This is also
known as equifinality.

Contact: Dr. Sofia Carozza, sofia.carozza[at]mrc.cbu.cam.ac.uk
%}
%%
function missclass = run_missclass(p)

% *** set repository directory ***
repo_dir = '/your/path/to/AdaptiveStochasticity';
% change directory to the data directory
cd(strcat(repo_dir,'/data/'));
% load global statistics
load('global_statistics_aal116_seed_400_matching.mat');
global_statistics = global_statistics_aal116_seed_400_matching; 
clear global_statistics_aal116_seed_400_matching
% set save directory
sdir = 'repo_dir/missclass_625';

missclass = zeros(1,625);

for i = 1:625
    classdata = [squeeze(global_statistics(:,p,:));
        squeeze(global_statistics(:,i,:))];
    class = [repelem('A',625) repelem('B',625)]';
    classifier = fitcsvm(classdata,class);
    cv = crossval(classifier);
    missclass(i) = kfoldLoss(cv);
    fprintf('completed %g comparison(s) for simulation %g...\n',i,p)
end

clear classdata class classifier cv

% change directory
cd(sdir);
% save file
save(sprintf('missclass_%g.mat',p),'missclass','-v7.3');
end
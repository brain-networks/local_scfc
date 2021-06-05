clear all
close all
clc

% functions come from:
%   • Brain connectivity toolbox (https://sites.google.com/site/bctnet/)
%   • Caio Seguin's github (https://github.com/caioseguin/connectomics/)
%   • Custom code written for this paper specifically

% if you use this script, please cite:
%
%   Esfahlani, F. Z., Faskowitz, J., Slack, J., Misic, B., & Betzel, R. 
%   (2021). Local structure-function relationships in human brain networks 
%   across the human lifespan. bioRxiv.
%   (https://doi.org/10.1101/2021.05.23.445128).

%% set up path structure and load data

% add functions to path
addpath(genpath('fcn'));

% load group-representative sc/fc/euclidean distance data
load data/hcp_group_mtrx-400.mat

% define weighted and binary versions of sc
a = sc;
abin = +(sc > 0);
n = length(sc);

%% generate predictor matrices
% this section generates a series of fully-weighted and signed matrices to
% be used as predictors for fc.

% for path-based measures, we transform connection weights to costs via the
% following equation: cost = weight.^-gamma. the user needs to specify the
% range of gamma values. similarly, for flow graphs the user needs to
% specify markov times at which to evaluate the matrix.

% gamma values
gammavals = [0.25,0.5,1,2];

% markov times
tvals = [1,2.5,5,10];

%% binary predictors
PLbin = distance_bin(abin);                     % path length
Gbin = expm(abin);                              % communicability
Cosbin = 1 - squareform(pdist(abin,'cosine'));  % cosine distance
SIbin = search_information(abin,'inv',false);   % search info
PTbin = path_transitivity(abin,'inv');          % path transitivity
mfptbin = zscore(mean_first_passage_time(abin));% mean first passage time
MIbin = matching_ind_und(abin);                 % matching index

FGbin = zeros(n,n,length(tvals));               % flow graphs
for itime = 1:length(tvals)
    FGbin(:,:,itime) = fcn_flow_graph(abin,ones(n,1),tvals(itime));
end

%% weighted predictors
Gwei = communicability_wei(a);                  % communicability
Coswei = 1 - squareform(pdist(a,'cosine'));     % cosine distance
mfptwei = zscore(mean_first_passage_time(a));   % mean first passage time
MIwei = matching_ind_und(a);                    % matching index
PLwei = zeros(n,n,length(gammavals));           % path length
SIwei = PLwei;                                  % search info
PTwei = PLwei;                                  % path transitivity
for igamma = 1:length(gammavals)
    
    L = a.^(-gammavals(igamma));                % convert weight to cost
    PLwei(:,:,igamma) = distance_wei_floyd(L);
    
    L(isinf(L)) = 0;
    SIwei(:,:,igamma) = search_information(L,[],false);
    PTwei(:,:,igamma) = path_transitivity(L,[]);
    
end
FGwei = FGbin;                                  % flow graphs
for itime = 1:length(tvals)
    FGwei(:,:,itime) = fcn_flow_graph(a,ones(n,1),tvals(itime));
    
end

%% navigation-based predictors
nav_struct = navigate(a,d);
failed = nav_struct.failed_paths;
nav_struct.num_hops(failed == 1) = inf;
nav_struct.pl_MS(failed == 1) = inf;

Navnumhops = nav_struct.num_hops;               % number of hops
NavplMS = nav_struct.pl_MS;                     % total distance

%% aggregate predictors and construct local multi-linear models
% concatenate all predictors into single array
mats = cat(3,d,PLbin,PLwei,Gwei,Gbin,Coswei,Cosbin,SIbin,SIwei,PTbin,PTwei,MIwei,MIbin,Navnumhops,NavplMS,mfptwei,mfptbin,FGbin,FGwei);

% preallocate array for storing correlation coefficients
r = zeros(n,size(mats,3));
for p = 1:size(mats,3)
    for i = 1:n
        
        % column of fc (from one region)
        y = fc(:,i);
        
        % column from predictor p
        x = mats(:,i,p);
        
        % remove node i
        y(i) = [];
        x(i,:) = [];
        
        % ignore nans/infs
        mask = ~isnan(x) & ~isinf(x) & ~isnan(y) & ~isinf(y);
        x = x(mask);
        y = y(mask);
        
        % z-score predictor
        x = x - nanmean(x);
        x = x/nanstd(x);
        
        % z-score fc
        y = y - nanmean(y);
        y = y/nanstd(y);
        
        % calculate correlation
        r(i,p) = corr(x,y);
        
        % note: could replace call to 'corr' with 'regress'
        %[beta,~,~,~,stats] = regress(y,[ones(length(y),1),x]);
        %s(i,p) = sign(beta(2))*sqrt(stats(1));
        
    end
end

% transform correlation coefficients into variance explained
rsq = r.^2;
%% A more statistically driven approach to sensitivity tests
% derives distributions for all parameters on which our model simulations
% will depend on: W, L, Zg, Zs, Ta, Sa, Da, Qa
run load_local_paths.m % sets data_path, import_path, collation_path, model_path, and project_path
addpath(genpath(import_path))
addpath(genpath(model_path))
addpath(genpath(collation_path))

outs_path = [data_path,'/greenland/FjordMIX/boxmodel/mcmc/']; % where the model output files will be saved
figs_path = [project_path,'/figs/mcmc/'];                     % where the figures and animations will be saved

[~,fjords_processed,~,~] = compile_datasets(data_path);

%% Gets distribution of fjord geometry parameters
verbose.print = 0;
verbose.plot  = 1;
fjord_stats = print_fjord_statistics(fjords_processed,verbose);

W_space  = fjord_stats.W.pd.random(10000,1);
L_space  = fjord_stats.L.pd.random(10000,1);
Zg_space = fjord_stats.Zg.pd.random(10000,1);
Zs_space = fjord_stats.Zs.pd.random(10000,1);
% figure; 
% subplot(2,2,1); histogram(L_space);
% subplot(2,2,2); histogram(W_space);
% subplot(2,2,3); histogram(Zg_space);
% subplot(2,2,4); histogram(Zs_space);

%% Sets up the ocean forcing

% Detrend T and S
% Compute T and S climatologies (Tclim, Sclim)
% Compute anomalies (Ta,Sa) = (T,S) - (Tclim,Sclim)
% Compute PDF of anomalies

%% Sets up the glacier forcing

% Subglacial discharge (same procedure as for T and S)

% Solid-ice discharge (same procedure as for T and S)

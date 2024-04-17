run load_local_paths.m % sets data_path, import_path, collation_path, model_path, and project_path
addpath(genpath(import_path))
addpath(genpath(model_path))
addpath(genpath(collation_path))
addpath(genpath(uqlab_path))
addpath(genpath('./'))

inputs_path = [data_path,'/greenland/Cowton2023GRL/'];                % where the input data is stored
outs_path   = [data_path,'/greenland/FjordMIX/boxmodel/cowton2023/']; % where the model output files will be saved
figs_path   = [project_path,'/figs/cowton2023/'];                     % where the figures and animations will be saved

fun = @(s) all(structfun(@isempty,s)); % tiny function to get rid of empty entries in array
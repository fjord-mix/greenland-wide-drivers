run load_local_paths.m % sets data_path, import_path, collation_path, model_path, and project_path
addpath(genpath(import_path))
addpath(genpath(model_path))
addpath(genpath(collation_path))
addpath(genpath(uqlab_path))
addpath(genpath('./'))

inputs_path = [data_path,'/greenland_common/Cowton2023GRL/'];                % where the input data is stored
outs_path   = '/Users/mmeb1/OneDrive - University of St Andrews/data_common/greenland/FjordMIX/rpm/'; % where the model output files will be saved
figs_path   = [project_path,'/figs/GRL_wide/'];                     % where the figures and animations will be saved

fun = @(s) all(structfun(@isempty,s)); % tiny function to get rid of empty entries in array
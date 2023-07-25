% RUN_BOXMODEL Driver to run the box model simulation.

%% Examples
% 
% 1. SGD-driven, constant shelf temp, stratified shelf salinity, vertical
% mixing.

% Explore parameter space
experiment = 'sgd';
parameter_to_vary = 'Qv0';
parameter_space = linspace(350, 1000, 50);
t = linspace(1,100,1);
if exist('fjord_par_outputs','Var'), clear fjord_par_outputs; end
fjord_par_outputs(length(parameter_space)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"s",[]);

%% Run the model (output is NOT saved automatically).
% Get run parameters

% parpool(2)
poolobj = gcp('nocreate');
addAttachedFiles(poolobj,regexp(genpath(model_path), ':', 'split'));
parfor INDEX = 1:length(parameter_space)    
    p_def = get_model_default_parameters;
    fjord_par_outputs(INDEX).m.name = sprintf('Iteration_%d\n', INDEX);
    fjord_par_outputs(INDEX).t = linspace(1,100,1);
    fjord_par_outputs(INDEX).p = mod_run_param(parameter_to_vary, parameter_space(INDEX), p_def);
    [fjord_par_outputs(INDEX).s,fjord_par_outputs(INDEX).f] = boxmodel(fjord_par_outputs(INDEX).p, t);
    disp(fjord_par_outputs(INDEX).m.name)
end
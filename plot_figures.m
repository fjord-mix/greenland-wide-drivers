function plot_figures(datasets, fjord, fjords_map, ctd_data)
% Model validation

% Loading model output
% file_fjord_run  = strsplit(input_file_name,'/');
% fname_fjord_run = strsplit(file_fjord_run{end},'.');
% fjord_title     = strsplit(fname_fjord_run{1},'_');
% fjord           = load(input_file_name).fjord_model_out;
fjord_title   = fjord.m.name;
model_runtime = fjord.s.t(1:size(fjord.s.H,2));
if nargin < 4
    ctd_data = datasets.obs.ctd_data.all;
end
if length(model_runtime) < 2
    disp('Cannot plot data, as there is only one recorded model time step!')
    return
end

%% Plotting figures

% Time series of T and S per layer of the desired fjord
[~] = plot_TS_boxmodel(datasets,fjords_map,fjord_title,fjord,model_runtime,ctd_data); % TODO: check if it still work
% print(['/Users/mmeb1/FjordMIX/presentations/2305_StAndrews/figs/out_ts_',fjord_title{2},'_shelf_only.png'],'-dpng','-r300');

% Time series of layer depths of the desired fjord
[~] = plot_H_boxmodel(fjord_title,fjord,model_runtime);
% print(['/Users/mmeb1/FjordMIX/presentations/2305_StAndrews/figs/out_h_',fjord_title{2},'_shelf_only'],'-dpng','-r300');


end % end function
n_years = 5;
res_box_yr     = cell([n_years,1]);
ensemble_yr    = cell(size(res_box_yr));
fjord_model_yr = cell(size(res_box_yr));
fjord_IDs      = char(65:1:65+length(fjord_names)-1);

for i_yr_load=1:n_years
    file_in = [outs_path,'rpm_NWfjords_n',num2str(n_runs),'_',num2str(2015+i_yr_load),'_',num2str(60),'layers_dt',num2str(dt_in_h),'h'];
    load(file_in);

    [~,res_box_yr{i_yr_load}] = postprocess_ensemble(fjord_model,ensemble,tgt_days);
    ensemble_yr{i_yr_load}    = ensemble;
    fjord_model_yr{i_yr_load} = fjord_model;
    fprintf('Postprocessing %d done.\n',2015+i_yr_load)
end
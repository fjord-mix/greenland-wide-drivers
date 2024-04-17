[datasets,fjords_compilation,fjords_map,~,glaciers_compilation] = compile_datasets(data_path);
datasets.opts.time_start = time_start;
datasets.opts.time_end   = time_end;
datasets.opts.dt         = 30; % time step (in days) for compiling the forcings (basically just a dummy here)
fjords_processed(size(fjords_compilation)) = struct("p",[],"a",[],"f",[],"t",[],"m",[]);
for i=1:length(fjords_compilation)
    fjords_processed(i) = prepare_boxmodel_input(datasets,fjords_compilation(i)); % yes, the code will complain, but we are not using the ocean data from it
end

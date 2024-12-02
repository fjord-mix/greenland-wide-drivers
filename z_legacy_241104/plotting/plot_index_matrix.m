function hf = plot_index_matrix(ax,analysisStruct,var,IOpts,regions,analysisString,cmap_n)

n_regions = length(regions);
n_inputs  = length(IOpts{1}.Marginals);
matrix_delta = NaN([n_inputs, n_regions]);

for i_reg=1:n_regions
    for i_input=1:n_inputs
        borgResults_ohc  = analysisStruct{i_reg}.Results;
        matrix_delta(i_input,i_reg) = borgResults_ohc.(var)(i_input);
    end
end

if isempty(ax)
    hf = figure('Name',analysisString,'position',[40 40 600 400]);
    ax = gca;
end
imagesc(ax,matrix_delta); hc1 = colorbar;

ylabel(hc1,analysisString,'FontSize',14,'FontName','DejaVu Sans')
set(ax,'YTick', 1:n_inputs, 'YTickLabel', borgResults_ohc.VariableNames,...
        'XTick', 1:n_regions,'XTickLabel',regions,...
        'FontSize',14,'FontName','DejaVu Sans');

if nargin < 7, cmap_n = cbrewer('seq','Reds',7,'linear'); end
colormap(ax,cmap_n)

end
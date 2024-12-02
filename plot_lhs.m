function hf = plot_lhs(X,param_names,param_units,log_axes)

if nargin < 4, log_axes=0; end

n_params = size(X,2);
if length(param_names) ~= n_params || length(param_units) ~= n_params
    error("Parameter names and units must agree with input matrix dimension")
end

combinations = nchoosek(1:n_params,2);

hf = figure('Position',[40 40 600 600]);
tiledlayout('flow');
for i_plt=1:size(combinations,1)
    nexttile; 
    hold on; box on; grid on;
    scatter(X(:,combinations(i_plt,1)),X(:,combinations(i_plt,2)),20,'k','filled')
    xlabel(sprintf('%s (%s)',param_names{combinations(i_plt,1)},param_units{combinations(i_plt,1)}))
    ylabel(sprintf('%s (%s)',param_names{combinations(i_plt,2)},param_units{combinations(i_plt,2)}))
    set(gca,'fontsize',14)
    if log_axes
        if (max(X(:,combinations(i_plt,1))) - min(X(:,combinations(i_plt,1))) > 1e3) || ...
            max(X(:,combinations(i_plt,1))) - min(X(:,combinations(i_plt,1))) < 1e-3
                set(gca,'XScale','log')
        end
        if (max(X(:,combinations(i_plt,2))) - min(X(:,combinations(i_plt,2))) > 1e3) || ...
            max(X(:,combinations(i_plt,2))) - min(X(:,combinations(i_plt,2))) < 1e-3
                set(gca,'YScale','log')
        end
    end
end


end
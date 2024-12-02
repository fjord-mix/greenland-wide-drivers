function hf = plot_inputs_summary(IOpts,probs)
% Plots all distributions to be used as input for the emulator/box model
% training set

hf = figure("Name",'Probability functions','Position',[50 50 1200 700]);
t=tiledlayout('flow');

for i=1:length(probs)
    nexttile; box on;
    opts = IOpts.Marginals(i);
    xi = linspace(0.95.*opts.Bounds(1),1.05.*opts.Bounds(2),100);
    prob = pdf(probs(i),xi);
    plot(xi,prob,'-k','linewidth',1.5)
    xlabel(opts.Name);
    xlim([0.95.*opts.Bounds(1),1.05.*opts.Bounds(2)])
    ylim([min(prob), 1.1.*max(prob)])
    set(gca,'FontSize',14)
end

ylabel(t,'Probability function','fontsize',14)
t.TileSpacing='compact';
t.Padding='compact';
end
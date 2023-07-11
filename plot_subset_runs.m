function plot_subset_runs(files,set_name,cmap_name)

f = dir(files);
n_files = length(f);
if nargin > 2
    line_colors=cbrewer('qual',cmap_name,n_files);
else
    line_colors=cbrewer('qual','Set1',n_files);
end


m=4; n=3; % 4 layers, H, T, and S
figure('Position',[20 20 900 900],'Name',set_name); 

labels={};
handles = [];
for i_file=1:n_files
    fjord_run=load([f(i_file).folder,'/',f(i_file).name]).fjord_run;
    labels{end+1} = f(i_file).name;
    model_runtime = fjord_run.s.t(1:size(fjord_run.s.H,2));
    runtime_axis = fjord_run.m.time_axis;
    t0 = convertTo(runtime_axis(1),'datenum');
    taxis = NaT([size(fjord_run.s.H,2),1]);
    for i_time=1:length(taxis)
        taxis(i_time) = datetime(t0+model_runtime(i_time),'ConvertFrom','datenum');
    end
    n_steps = length(taxis);

    n_layers=size(fjord_run.s.H,1);

    %% Thickness
    for i_layer=1:n_layers
        subplot(m,n,1+(i_layer-1)*3); hold on; box on
        if i_layer==1, title(sprintf('Layer thickness (m)')); end
        hp = plot(taxis,real(fjord_run.s.H(i_layer,1:n_steps)),'linewidth',1.5,'color',line_colors(i_file,:)); 
        text(0.02,0.9,sprintf('Layer %d',i_layer),'Units','normalized')    
        if i_layer==4          
            yticks(mean(fjord_run.s.H(i_layer,1:n_steps)))
            yticklabels(sprintf('%0.2f',mean(fjord_run.s.H(i_layer,1:n_steps))));
        end
    end
    handles = [handles; hp];
    

    %% Temperature
    for i_layer=1:n_layers
        subplot(m,n,2+(i_layer-1)*3); hold on; box on
        if i_layer==1, title(sprintf('Temperature ($^o$C)','interpreter','latex')); end
        plot(taxis,real(fjord_run.s.T(i_layer,1:n_steps)),'linewidth',1.5,'color',line_colors(i_file,:)); 
        % text(0.02,0.9,sprintf('Layer %d',i_layer),'Units','normalized')
    end
    
    %% Salinity
    for i_layer=1:n_layers
        subplot(m,n,3+(i_layer-1)*3); hold on; box on
        if i_layer==1, title(sprintf('Salinity (-)')); end
        plot(taxis,real(fjord_run.s.S(i_layer,1:n_steps)),'linewidth',1.5,'color',line_colors(i_file,:));     
        % text(0.02,0.9,sprintf('Layer %d',i_layer),'Units','normalized')
    end
    
end
% legend('plume','shelf','mixing','artificial','icebergs','Location','east');
legend(handles,labels,'interpreter','none','Location','east')


end
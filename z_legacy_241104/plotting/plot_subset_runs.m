function plot_subset_runs(files,set_name,cmap_name)

isdays='';
f = dir(files);
n_files = length(f);
if nargin > 2
    line_colors=cbrewer('qual',cmap_name,n_files);
else
    line_colors=cbrewer('qual','Set1',n_files);
end
labels={};
handles = [];

fjord_run=load([f(1).folder,'/',f(1).name]).fjord_run;
if isfield(fjord_run,'o')
    figure('Position',[20 20 800 600],'Name',set_name); 

    for i_file=1:n_files
        fjord_run=load([f(i_file).folder,'/',f(i_file).name]).fjord_run;        
        if ~isfield(fjord_run,'o'), continue; end % will skip if no postprocessed structure exists
        model_runtime = fjord_run.s.t(1:size(fjord_run.s.H,2));

        if isfield(fjord_run.m,'time_axis')
            runtime_axis = fjord_run.m.time_axis;
            t0 = convertTo(runtime_axis(1),'datenum');
            taxis = NaT([size(fjord_run.s.H,2),1]);
            for i_time=1:length(taxis)
                taxis(i_time) = datetime(t0+model_runtime(i_time),'ConvertFrom','datenum');
            end
        else
            taxis=fjord_run.s.t;
            isdays = ' (days)';
        end
        

        subplot(2,1,1); hold on; box on
        hp = plot(taxis,sum(fjord_run.o.hc),'linewidth',1.5,'color',line_colors(i_file,:)); 
        ylabel('Heat content (J)','interpreter','latex'); 

        subplot(2,1,2); hold on; box on
        plot(taxis,sum(fjord_run.o.sc),'linewidth',1.5,'color',line_colors(i_file,:)); 
        ylabel('Salt content (kg)','interpreter','latex');
        xlabel(['Model time',isdays])

        labels{end+1} = f(i_file).name;
        handles       = [handles; hp];
    end
else
    m=4; n=3; % 4 layers, H, T, and S
    figure('Position',[20 20 900 900],'Name',set_name); 
    
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
end
legend(handles,labels,'interpreter','none','Location','southeast')

end
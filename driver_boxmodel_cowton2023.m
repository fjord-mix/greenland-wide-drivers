%% Driver file for simulating the fjords presented in Cowton et al. (2023, GRL)

%% Configuring paths
run load_local_paths.m % sets data_path, import_path, collation_path, model_path, and project_path
addpath(genpath(import_path))
addpath(genpath(model_path))
addpath(genpath(collation_path))
addpath(genpath(uqlab_path))
addpath(genpath('./'))

inputs_path = [data_path,'/greenland/Cowton2023GRL/'];                % where the input data is stored
outs_path   = [data_path,'/greenland/FjordMIX/boxmodel/cowton2023/']; % where the model output files will be saved
figs_path   = [project_path,'/figs/cowton2023/'];                     % where the figures and animations will be saved


fun = @(s) all(structfun(@isempty,s)); % tiny function to get rid of empty entries in array
%% Compiling data for all fjords around Greenland (requires data-compilation repository)
% [datasets,fjords_compilation,fjords_map,~,glaciers_compilation] = compile_datasets(data_path);
% datasets.opts.time_start = time_start;
% datasets.opts.time_end   = time_end;
% datasets.opts.dt         = 30; % time step (in days) for compiling the forcings (basically just a dummy here)
% fjords_processed(size(fjords_compilation)) = struct("p",[],"a",[],"f",[],"t",[],"m",[]);
% for i=1:length(fjords_compilation)
%     fjords_processed(i) = prepare_boxmodel_input(datasets,fjords_compilation(i));
% end
% 
% % These are the IDs of the corresponding fjords above in the "fjords_processed" data structure
% fjord_ids = [4,9,17,20,22,23,24,25,28,29,30,31,24,37];


%% Initialise all needed variables

i_yr = 4; % we want only 2019 data for now, because we only have geometry data for glaciers with data for 2019
time_start = datetime(2019,01,01);
time_end   = datetime(2019,12,31);
input_dt   = 30;
model_dt   = 0.1;

letters = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n'};

fjord_names = {'Qeqertaarsuusarsuup','Nuussuup','Tuttulikassaap','Nunatakassaap','Illulip',...
               'Kakiffaat','Naajarsuit','Sermeq',...
               'Salliarsutsip','Umiammakku','Kangilliup','Kangerlussuup','SermeqSilarleq','SermeqKujalleq'} ;


%% Reading the data from Cowton et al. (2023)

% modelled subglacial discharge
qsg_glaciers = load([inputs_path,'/runoff/runoff.txt']);
time_glaciers = datetime(load([inputs_path,'runoff/runoff_tvec.txt']));
tgt_period = time_glaciers > time_start & time_glaciers < time_end;
qsg_all = qsg_glaciers(tgt_period,:);
taxis_qsg = time_glaciers(tgt_period);
clear qsg_glaciers time_glaciers tgt_period % bit of tidying up

% CTD profiles
ctd_data   = load([inputs_path,'/fjords/NW_fjord_TSz']);
meta_table = readtable([inputs_path,'/fjords/Fjords_table.xlsx'],'Range','A3:I17'); % I had to unmerge some sells from the original table, otherwise Matlab would not get the header names (why were they merged in the first place?!)



%% Preparing the model runs
n_fjords = size(meta_table,1);

if exist('fjord_model',"var"),       clear fjord_model; end
fjord_model(n_fjords) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[]);
for i_fjord=1:n_fjords
    if ~isnan(meta_table.Length_km_(i_fjord)) % we only want fjords for which we have all geometry data

        % time axis for the run
        t = 0:model_dt:3650; % 10 years (we want to repeat the yearly data 10 times)

        % compile metadata and default parameters
        m.ID   = meta_table.Var1(i_fjord);
        m.name = meta_table.GreenlandicName(i_fjord);
        m.lon  = meta_table.Longitude__E_(i_fjord);
        m.lat  = meta_table.Latitude__N_(i_fjord);
        [p,~]  = get_model_default_parameters();
        
        % obtain CTD data
        ctd.temp  = ctd_data.Tdata{i_yr,i_fjord};
        ctd.salt  = ctd_data.Sdata{i_yr,i_fjord};
        ctd.depth = ctd_data.zdata{i_yr,i_fjord};
        

        % fjord geometry
        p.L         = 1e3.*meta_table.Length_km_(i_fjord);
        p.W         = 1e3.*meta_table.Width_km_(i_fjord);
        p.silldepth = -1.* meta_table.SillDepth_m_(i_fjord);
        p.zgl       = -1.* meta_table.GroundingLineDepth_m_(i_fjord);
        if p.silldepth < p.zgl % if the GL sits above the sill
            p.H = abs(p.silldepth-p.Hmin); % we consider the fjord depth to be just Hmin below the sill, 
                                      % so we can avoid fiddling with the number of layers
        else
            p.H = abs(p.zgl); % if the GL is below the sill, then we treat it as the same as the fjord depth
        end

        % ocean forcing
        ts = ctd.temp.s;
        ss = ctd.salt.s;
        zs = -ctd.depth.s;
        if max(abs(zs)) < p.H % if our profile does not extend to the bottom of the fjord
            h_missing = p.H - max(abs(zs));
            depth_range_missing = max(abs(zs))+1:1:p.H;
            new_zs = [zs -depth_range_missing];
            ts = interp1(zs,ts,new_zs,'nearest','extrap');
            ss = interp1(zs,ss,new_zs,'nearest','extrap');
            zs = new_zs;
            clear new_zs h_missing depth_range_missing % tidying up
        end
        f.zs = zs;
        f.Ts = repmat(ts,length(t),1);
        f.Ss = repmat(ss,length(t),1);


        % initial conditions
        H_clim = double(get_fjord_boxes_from_density(f.Ts',f.Ss',f.zs,p)); % negative values, top to bottom
        [T_clim,S_clim] = bin_ocean_profiles(mean(f.Ts,1),mean(f.Ss,1),f.zs,H_clim,p);
        a.H0 = H_clim;
        a.T0 = T_clim;
        a.S0 = S_clim;

        % flipping the ocean forcing like the model wants it
        f.zs = flip(f.zs);
        f.Ts = flip(f.Ts,2)'; 
        f.Ss = flip(f.Ss,2)';
        
        % Sanity-check plot to make sure the boxes and binned T and S make sense
        % figure;
        % ints=[0,-cumsum(a.H0)];
        % subplot(1,2,1); hold on
        % plot(ts,zs);
        % scatter(a.T0,(ints(1:end-1)+ints(2:end))/2,'filled');
        % for j=1:size(ints,2)
        %     yline(ints(j));
        % end
        % yline(p.zgl,'--k','linewidth',0.5)
        % yline(p.silldepth,'--k','linewidth',0.5)
        % ylim([-p.H 0])
        % subplot(1,2,2); hold on
        % plot(ss,zs);
        % scatter(a.S0,(ints(1:end-1)+ints(2:end))/2,'filled');
        % for j=1:size(ints,2)
        %     yline(ints(j));
        % end
        % yline(p.zgl,'--k','linewidth',0.5)
        % yline(p.silldepth,'--k','linewidth',0.5)
        % ylim([-p.H 0])
        
        % obtain Qsg
        % interpolate from daily Qsg data to the model time axes
        taxis_qsg_julian = juliandate(taxis_qsg);
        taxis_qsg_num = taxis_qsg_julian-taxis_qsg_julian(1);
        taxis_qsg_num_dtmodel = taxis_qsg_num(1):t(2)-t(1):taxis_qsg_num(end);
        qsg_dtmodel = interp1(taxis_qsg_num,qsg_all(:,i_fjord),taxis_qsg_num_dtmodel);
        f.Qsg = repmat(qsg_dtmodel,1,11);
        f.Qsg = f.Qsg(1:length(t));
        clear taxis_qsg_julian taxis_qsg_num taxis_qsg_num_dtmodel qsg_dtmodel % bit of tidying up

        % no icebergs for now
        p.M0=0; 
        f.zi = -800:10:0;
        f.zi = f.zi';
        f.D  = zeros(1,length(t));
        f.xi = (p.nu0/sum(a.H0))*exp(p.nu0*f.zi/sum(a.H0))/(1-exp(-p.nu0));
        a.I0 = 0*f.zi;
        
        % obtain iceberg data - for later
        % runtime_axis =
        % datasets.opts.time_start:model_dt:datasets.opts.time_end;
        % [f.zi,f.D,f.xi,a.I0] = get_total_iceberg_discharge(datasets,fjords_compilation(fjord_ids(i_fjord)).glaciers,runtime_axis',p,a);
        % we need to repeat the 1-year discharge for 10 years like with Qsg

        fjord_model(i_fjord).t = t;
        fjord_model(i_fjord).m = m;
        fjord_model(i_fjord).p = p;
        fjord_model(i_fjord).a = a;
        fjord_model(i_fjord).f = f;
        fjord_model(i_fjord).c.tf = ctd.temp.f;
        fjord_model(i_fjord).c.sf = ctd.salt.f;
        fjord_model(i_fjord).c.zf = ctd.depth.f;
        fjord_model(i_fjord).c.ts = ctd.temp.s;
        fjord_model(i_fjord).c.ss = ctd.salt.s;
        fjord_model(i_fjord).c.zs = ctd.depth.s;
        
    end
end

idx = arrayfun(fun,fjord_model);
fjord_model(idx)=[]; % remove the empty elements

%% Run the model

n_fjord_runs = length(fjord_model);
range_C0 = [1e3,5e3,1e4,5e4,1e5];
range_P0 = 5:5:30;
n_C0 = length(range_C0);
n_P0 = length(range_P0);
if exist('ensemble',"var"),clear ensemble; end
% ensemble(n_fjord_runs, length(range_C0), length(range_P0)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[],"s",[]);
ensemble(n_fjord_runs, length(range_C0), length(range_P0)) = struct("p",[],"t",[],"m",[],"s",[]);
run_counter=0;
for i_fjord=1:n_fjord_runs
    for i_C0=1:n_C0
        for i_P0=1:n_P0
            run_counter = run_counter+1;
            cur_fjord = fjord_model(i_fjord);
            cur_fjord.p.C0 = range_C0(i_C0);
            cur_fjord.p.P0 = range_P0(i_P0);
            % cur_fjord.p.plot_runtime=1;
            tic
            try
                [cur_fjord.s,cur_fjord.f] = boxmodel(cur_fjord.p, cur_fjord.t, cur_fjord.f, cur_fjord.a);
                fprintf('run %d complete\n',run_counter)
                ensemble(i_fjord,i_C0,i_P0).p = cur_fjord.p;
                ensemble(i_fjord,i_C0,i_P0).t = cur_fjord.t;
                ensemble(i_fjord,i_C0,i_P0).m = cur_fjord.m;
                ensemble(i_fjord,i_C0,i_P0).s.Tfinal = mean(cur_fjord.s.T(:,end-365:end),2);
                ensemble(i_fjord,i_C0,i_P0).s.Sfinal = mean(cur_fjord.s.S(:,end-365:end),2);
                ensemble(i_fjord,i_C0,i_P0).s.Hfinal = mean(cur_fjord.s.H(:,end-365:end),2);
            catch ME
                fprintf('run %d failed: %s\n',run_counter,ME.message)
            end
            toc
        end
    end
end
save([outs_path,'boxmodel_runs_n',num2str(n_C0*n_P0)],'-v7.3','ensemble','fjord_model');
%% post-processing to make things 1:1 comparable

if exist('res_obs',"var"),       clear res_obs; end
res_obs(size(ensemble)) = struct("tf",[],"sf",[],"zf",[]);
if exist('res_box',"var"),       clear res_box; end
res_box(size(ensemble)) = struct("tf",[],"sf",[],"zf",[],"ID",[],"name",[]);
for i=1:size(ensemble,1)
    zf_obs = fjord_model(i).c.zf;
    tf_obs = fjord_model(i).c.tf;
    sf_obs = fjord_model(i).c.sf;
    n_completed = 0; % completed runs for that fjord

    % n_layers = ensemble(i,1,1).p.N+ensemble(i,1,1).p.sill;
    tf_box_comp = NaN([length(zf_obs),n_fjord_runs,length(range_C0),length(range_P0)]);
    sf_box_comp = NaN(size(tf_box_comp));
    for j=1:size(ensemble,2)
        for k=1:size(ensemble,3)
            if ~isempty(ensemble(i,j,k).s)
                ints=[0; cumsum(ensemble(i,j,k).s.Hfinal)];
                z_box = (ints(1:end-1)+ints(2:end))/2;
                tf_box = ensemble(i,j,k).s.Tfinal; 
                sf_box = ensemble(i,j,k).s.Sfinal; 
                
                tf_box_comp(:,i,j,k) = interp1(z_box,tf_box,zf_obs,'nearest','extrap');
                sf_box_comp(:,i,j,k) = interp1(z_box,sf_box,zf_obs,'nearest','extrap');
                n_completed = n_completed+1;
            end
        end
    end
    res_obs(i).zf = zf_obs;
    res_obs(i).tf = tf_obs;
    res_obs(i).sf = sf_obs;
    
    res_obs(i).zs = fjord_model(i).c.zs;
    res_obs(i).ts = fjord_model(i).c.ts;
    res_obs(i).ss = fjord_model(i).c.ss;

    res_box(i).zf = z_box;
    res_box(i).tf = mean(tf_box_comp(:,i,:,:),[3, 4],'omitnan');
    res_box(i).tfmin = min(tf_box_comp(:,i,:,:),[],[3, 4],'omitnan');
    res_box(i).tfmax = max(tf_box_comp(:,i,:,:),[],[3, 4],'omitnan');
    res_box(i).sf = mean(sf_box_comp(:,i,:,:),[3, 4],'omitnan');
    res_box(i).id = fjord_model(i).m.ID{1};
    res_box(i).name = fjord_model(i).m.name{1};
    res_box(i).n = n_completed;
end

% idx = arrayfun(fun,ensemble);
% ensemble(idx)=[]; % remove the empty elements

%% Plotting results
lcolor = lines(3);
figure('Name','Temperature comparison','Position',[40 40 1200 800]);
for i=1:n_fjord_runs
    subplot(2,5,i); hold on; box on; grid on
    text(0.02,1.05,sprintf("(%s) %s",res_box(i).id,res_box(i).name),'units','normalized','fontsize',14)
    text(0.02,0.05,sprintf("n=%d",res_box(i).n),'Units','normalized','FontSize',14)

    % figure; hold on
    % y2 = [res_obs(i).zf; flip(res_obs(i).zf)]';
    % inBetween = [res_box(i).tfmin, flip(res_box(i).tfmax)];
    % hp = fill(flip(inBetween,1), flip(-y2,2), lcolor(3,:),'edgecolor','none','facealpha',0.2);

    hs = plot(res_obs(i).ts,-res_obs(i).zs,'linewidth',1.5,'color',lcolor(1,:));
    hf = plot(res_obs(i).tf,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(2,:));
    hm = plot(res_box(i).tf,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:));
    plot(res_box(i).tfmin,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
    plot(res_box(i).tfmax,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
    if ismember(i,[1,6]),ylabel('Depth (m)'); end
    if i > 5, xlabel('Temperature (^oC)'); end
    set(gca,'fontsize',14)
    xlim([-1 6])
end
legend([hs, hf, hm], {"Shelf","Fjord","Model"},'fontsize',14,'Location','Southeast')

exportgraphics(gcf,[figs_path,'temp_profiles_2019_n',num2str(n_C0*n_P0),'.png'],'Resolution',300)
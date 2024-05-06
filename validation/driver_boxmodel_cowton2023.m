%% Driver file for simulating the fjords presented in Cowton et al. (2023, GRL)

%% Configuring paths
run setup_paths
fun = @(s) all(structfun(@isempty,s)); % tiny function to get rid of empty entries in array

%% Initialise all needed variables

i_yr = 4; % {1,2,3,4} for {2016,2017,2018,2019}
time_start = datetime(2015+i_yr,01,01);
time_end   = datetime(2015+i_yr,12,31);
input_dt   = 30;
model_dt   = 0.1;
n_years    = 4; % how many years we want to run
tgt_days   = [935,1010,1150]; % which days of the run we want to compare: year 3 peak-melt and post-melt, following pre-melt
name_days  = {'peak','post','winter'};

letters = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n'};

fjord_names = {'Qeqertaarsuusarsuup','Nuussuup','Tuttulikassaap','Nunatakassaap','Illulip',...
               'Kakiffaat','Naajarsuit','Sermeq',...
               'Salliarsutsip','Umiammakku','Kangilliup','Kangerlussuup','SermeqSilarleq','SermeqKujalleq'} ;

%% Compiling data 

run compile_process_fjords % Read compilation of all fjords around Greenland (requires data-compilation repository)

% These are the IDs of the corresponding fjords above in the "fjords_processed" data structure
fjord_ids = [4,9,17,20,22,23,24,25,28,29,30,31,24,37];

run load_cowton2023_data % Reading the data from Cowton et al. (2023)

tgt_period = time_glaciers > time_start & time_glaciers < time_end;
qsg_all = qsg_glaciers(tgt_period,:);
taxis_qsg = time_glaciers(tgt_period);
clear qsg_glaciers time_glaciers tgt_period % bit of tidying up

%% Preparing the model runs
n_fjords = size(meta_table,1);

if exist('fjord_model',"var"),       clear fjord_model; end
fjord_model(n_fjords) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[]);
for i_fjord=1:n_fjords
    if ~isnan(meta_table.Length_km_(i_fjord)) % we only want fjords for which we have all geometry data

        % time axis for the run
        t = 0:model_dt:365*n_years; % 4 years (we want to repeat the yearly data 4 times)

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
        
        % obtain Qsg
        % interpolate from daily Qsg data to the model time axes
        taxis_qsg_julian = juliandate(taxis_qsg);
        taxis_qsg_num = taxis_qsg_julian-taxis_qsg_julian(1);
        taxis_qsg_num_dtmodel = taxis_qsg_num(1):t(2)-t(1):taxis_qsg_num(end);
        qsg_dtmodel = interp1(taxis_qsg_num,qsg_all(:,i_fjord),taxis_qsg_num_dtmodel);
        f.Qsg = repmat(qsg_dtmodel,1,11);
        f.Qsg = f.Qsg(1:length(t));
        clear taxis_qsg_julian taxis_qsg_num taxis_qsg_num_dtmodel qsg_dtmodel % bit of tidying up

        % no icebergs
        % p.M0=0; 
        % f.zi = -800:10:0;
        % f.zi = f.zi';
        % f.D  = zeros(1,length(t));
        % f.xi = (p.nu0/abs(p.zgl))*exp(p.nu0*f.zi/abs(p.zgl))/(1-exp(-p.nu0));
        % a.I0 = 0*f.zi;
        
        % obtain iceberg data
        runtime_axis = datasets.opts.time_start:model_dt:datasets.opts.time_end;
        [f.zi,ice_d,f.xi,a.I0] = get_total_iceberg_discharge(datasets,fjords_compilation(fjord_ids(i_fjord)).glaciers,runtime_axis',p,a);
        % we need to repeat the 1-year discharge for 10 years like with Qsg
        f.D = repmat(ice_d,1,n_years+1);
        f.D = f.D(1:length(t));

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

        % Sanity-check plot to make sure the boxes and binned T and S make sense
        % % Sanity-check plot to make sure the boxes, T,S profiles, Qsg and D all make sense
        figure('Position',[100 100 800 500]);
        ints=[0,-cumsum(a.H0)];
        subplot(1,3,1); hold on; box on;
        plot(f.Ts(:,1),f.zs);
        scatter(a.T0,(ints(1:end-1)+ints(2:end))/2,'filled');
        for j=1:size(ints,2)
            yline(ints(j));
        end
        yline(p.zgl,'--k','linewidth',0.5)
        yline(p.silldepth,'--k','linewidth',0.5)
        ylim([-p.H 0])
        subplot(1,3,2); hold on; box on;
        plot(f.Ss(:,1),f.zs);
        scatter(a.S0,(ints(1:end-1)+ints(2:end))/2,'filled');
        for j=1:size(ints,2)
            yline(ints(j));
        end
        yline(p.zgl,'--k','linewidth',0.5)
        yline(p.silldepth,'--k','linewidth',0.5)
        ylim([-p.H 0])
        subplot(1,3,3); hold on; box on; grid on;
        plot(t,f.Qsg);
        plot(t,f.D);
        legend('Subglacial discharge','Solid-ice discharge');
        exportgraphics(gcf,[figs_path,'forcing_summary_c13_fjord_',letters{i_fjord},'_',num2str(2015+i_yr),'.png'],'Resolution',300)
        
        clear m p a f % bit of tidying up
    end
end

idx = arrayfun(fun,fjord_model);
fjord_model(idx)=[]; % remove the empty elements
disp('Model inputs processing done.')

%% Run the model

n_fjord_runs = length(fjord_model);
range_C0 = [1e3,5e3,1e4,5e4,1e5];
range_P0 = 5:5:30;
range_M0 = [5e-9,1e-8,2e-8,5e-8];
range_gamma = 0.4:0.05:0.6;
n_C0 = length(range_C0);
n_P0 = length(range_P0);
n_M0 = length(range_M0);
n_gamma = length(range_gamma);
n_combinations=n_gamma*n_M0*n_P0*n_C0;
if exist('ensemble',"var"),clear ensemble; end
% ensemble(n_fjord_runs, length(range_C0), length(range_P0)) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[],"s",[]);
ensemble(n_fjord_runs, length(range_C0), length(range_P0),length(range_M0),length(range_gamma)) = struct("p",[],"t",[],"m",[],"s",[]);
run_counter=0;
% i_M0=1; i_gamma=1;
% i_P0=1;
for i_fjord=1:n_fjord_runs
    tic
    zf_obs = fjord_model(i_fjord).c.zf;
    tf_obs = fjord_model(i_fjord).c.tf;
    sf_obs = fjord_model(i_fjord).c.sf;
    
    for i_C0=1:n_C0
    for i_P0=1:n_P0
    for i_M0=1:n_M0
    for i_gamma=1:n_gamma
        run_counter = run_counter+1;
        cur_fjord = fjord_model(i_fjord);
        cur_fjord.p.C0 = range_C0(i_C0);
        cur_fjord.p.P0 = range_P0(i_P0);
        cur_fjord.p.M0 = range_M0(i_M0);
        cur_fjord.p.gamma = range_gamma(i_gamma);
        % cur_fjord.p.plot_runtime=1;
        try
            [cur_fjord.s,cur_fjord.f] = boxmodel(cur_fjord.p, cur_fjord.t, cur_fjord.f, cur_fjord.a);
            fprintf('run %d complete. ',run_counter)
    
            ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).p = cur_fjord.p;
            ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).t = cur_fjord.t;
            ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).m = cur_fjord.m;
            Tfinal=NaN([cur_fjord.p.N+cur_fjord.p.sill,length(tgt_days)]);
            Sfinal=NaN(size(Tfinal));
            Hfinal=NaN(size(Tfinal));
            for i_day=1:length(tgt_days)
                % 10-day avg centered at the target day
                tgt_day = tgt_days(i_day);
                Tfinal(:,i_day) = mean(cur_fjord.s.T(:,(tgt_day-5:tgt_day+5)/model_dt),2); 
                Sfinal(:,i_day) = mean(cur_fjord.s.S(:,(tgt_day-5:tgt_day+5)/model_dt),2);
                Hfinal(:,i_day) = mean(cur_fjord.s.H(:,(tgt_day-5:tgt_day+5)/model_dt),2);
            end
            ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tfinal = Tfinal;
            ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Sfinal = Sfinal;
            ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Hfinal = Hfinal;

            Tregular=NaN([length(zf_obs),length(cur_fjord.s.t)]);
            for k=1:length(cur_fjord.s.t)
                ints=[0; cumsum(cur_fjord.s.H(:,k))];
                z_box = (ints(1:end-1)+ints(2:end))/2;
                Tregular(:,k) = interp1(z_box,cur_fjord.s.T(:,k),zf_obs,'nearest','extrap');
            end
            range_upper = zf_obs < 50;
            range_inter = zf_obs > 50 & zf_obs < 250;
            range_lower = zf_obs > 250 & zf_obs < 500;
            
            ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tupper = mean(Tregular(range_upper,:),1,'omitnan');   % avg temp of the upper 50 m
            ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tinter = mean(Tregular(range_inter,:),1,'omitnan');  % avg temp of the 50 - 250 m layer
            ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tlower = mean(Tregular(range_lower,:),1,'omitnan'); % avg temp of the 250 - 500 m layer
            fprintf('Output interpolation complete. ')
        catch ME
            fprintf('run %d failed: %s. ',run_counter,ME.message)
        end
        fprintf('\n')
    end % for i_gamma
    end % for i_M0
    end % for i_P0
    end % for i_C0
fprintf('Done with fjord %d. ',i_fjord)
toc
fprintf('\n')
end % for i_fjord
save([outs_path,'boxmodel_runs_c13comp_n',num2str(n_combinations),'_',num2str(2015+i_yr)],'-v7.3','ensemble','fjord_model');
disp('Outputs saved.')
% load([outs_path,'boxmodel_runs_MITgcm_obsice_comp_n',num2str(n_combinations),'_',num2str(2015+i_yr)]);

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
    tf_box_comp = NaN([length(zf_obs),n_fjord_runs,n_C0,n_P0,n_M0,n_gamma,length(tgt_days)]);
    sf_box_comp = NaN(size(tf_box_comp));
    tupper_box_comp = NaN([length(fjord_model(i).t)-1,n_fjord_runs,n_C0,n_P0,n_M0,n_gamma]);
    tinter_box_comp = NaN(size(tupper_box_comp));
    tlower_box_comp = NaN(size(tupper_box_comp));
    for j=1:size(ensemble,2)
    for k=1:size(ensemble,3)
    for l=1:size(ensemble,4)
    for n=1:size(ensemble,5) % we skip 'm' because it is being used for the metadata structure elsewhere in the code
        if ~isempty(ensemble(i,j,k,l,n).s) & ~isnan(ensemble(i,j,k,l,n).s.Hfinal(1))
            for i_day=1:length(tgt_days) % we repeat the whole thing for each of our target days
                ints=[0; cumsum(ensemble(i,j,k,l,n).s.Hfinal(:,i_day))];
                z_box = (ints(1:end-1)+ints(2:end))/2;
                tf_box = ensemble(i,j,k,l,n).s.Tfinal(:,i_day); 
                sf_box = ensemble(i,j,k,l,n).s.Sfinal(:,i_day);
                
                tf_box_comp(:,i,j,k,l,n,i_day) = interp1(z_box,tf_box,zf_obs,'nearest','extrap');
                sf_box_comp(:,i,j,k,l,n,i_day) = interp1(z_box,sf_box,zf_obs,'nearest','extrap');
            end

            tupper_box_comp(:,i,j,k,l,n) = ensemble(i,j,k,l,n).s.Tupper;
            tinter_box_comp(:,i,j,k,l,n) = ensemble(i,j,k,l,n).s.Tinter;
            tlower_box_comp(:,i,j,k,l,n) = ensemble(i,j,k,l,n).s.Tlower;

            n_completed = n_completed+1;
        end
    end % n
    end % l
    end % k
    end % j
    res_obs(i).zf = zf_obs;
    res_obs(i).tf = tf_obs;
    res_obs(i).sf = sf_obs;
    
    res_obs(i).zs = fjord_model(i).c.zs;
    res_obs(i).ts = fjord_model(i).c.ts;
    res_obs(i).ss = fjord_model(i).c.ss;

    res_box(i).t  = t(2:end); % ensemble(i,1,1,1,1).t(2:end); % need to come up with a good way to avoid a crashed run with no time variable
    res_box(i).zf = z_box;
    res_box(i).tf = median(tf_box_comp(:,i,:,:,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).tfmin = min(tf_box_comp(:,i,:,:,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).tfmax = max(tf_box_comp(:,i,:,:,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).sf = median(sf_box_comp(:,i,:,:,:,:,:),[3,4,5,6],'omitnan');


    res_box(i).Tupper = median(tupper_box_comp(:,i,:,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).Tupper_min = min(tupper_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).Tupper_max = max(tupper_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');

    res_box(i).Tinter = median(tinter_box_comp(:,i,:,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).Tinter_min = min(tinter_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).Tinter_max = max(tinter_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');

    res_box(i).Tlower = median(tlower_box_comp(:,i,:,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).Tlower_min = min(tlower_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).Tlower_max = max(tlower_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');

    res_box(i).ID = fjord_model(i).m.ID{1};
    res_box(i).name = fjord_model(i).m.name{1};
    res_box(i).n = n_completed;
end
disp('Postprocessing done.')

%% Plotting results
lcolor = lines(3+length(tgt_days));
fsize=12;
hf_profiles = figure('Name','Temperature profiles','Position',[40 40 1200 400*length(fjord_model)/2]);
tiledlayout(length(fjord_model)/2,2);
hf_series   = figure('Name','Temperature evolution','Position',[40 40 1200 400*length(fjord_model)/2]);
tiledlayout(length(fjord_model)/2,2);
for i=1:n_fjord_runs

    figure(hf_profiles)
    nexttile; hold on; box on; grid on
    text(0.02,1.02,sprintf("(%s) %s (%.0f km long)",res_box(i).ID,res_box(i).name, fjord_model(i).p.L/1e3),'units','normalized','VerticalAlignment','bottom','fontsize',fsize)
    text(0.02,0.02,sprintf("n=%d",res_box(i).n),'Units','normalized','VerticalAlignment','bottom','FontSize',fsize)

    % figure; hold on
    % y2 = [res_obs(i).zf; flip(res_obs(i).zf)]';
    % inBetween = [res_box(i).tfmin, flip(res_box(i).tfmax)];
    % hp = fill(flip(inBetween,1), flip(-y2,2), lcolor(3,:),'edgecolor','none','facealpha',0.2);

    hs = plot(res_obs(i).ts,-res_obs(i).zs,'linewidth',1.5,'color',lcolor(1,:));
    hf = plot(res_obs(i).tf,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(2,:));
    hb=[];
    for i_day=1:length(tgt_days)
        hb_d = plot(res_box(i).tf(:,i_day),-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3+i_day,:));
        hb = [hb hb_d];
    end
    % plot(res_box(i).tfmin,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
    % plot(res_box(i).tfmax,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
    scatter(0,fjord_model(i).p.zgl,40,'v','filled','MarkerFaceColor','black')
    plot([0 0],[-fjord_model(i).p.H fjord_model(i).p.silldepth],'-k','linewidth',2)

    if mod(i,2) > 0, ylabel('Depth (m)'); end
    if i>n_fjord_runs-2
        xlabel('Temperature (^oC)');  
    else
        set(gca,'xticklabels',[])
    end
    set(gca,'fontsize',fsize)
    xlim([-2 6])
    ylim([-fjord_model(i).p.H 0])

    figure(hf_series)
    nexttile; hold on; box on; grid on
    text(0.02,0.99,sprintf("(%s) %s (%.0f km long; n=%d)",res_box(i).ID,res_box(i).name, fjord_model(i).p.L/1e3,res_box(i).n),'units','normalized','VerticalAlignment','top','fontsize',fsize)
    % text(0.02,0.02,sprintf("n=%d",res_box(i).n),'Units','normalized','VerticalAlignment','bottom','FontSize',fsize)
    if length(res_box(i).t) == length(res_box(i).Tupper)
        
        hu = plot(res_box(i).t,res_box(i).Tupper,'linewidth',1.5,'color',lcolor(1,:));
        hi = plot(res_box(i).t,res_box(i).Tinter,'linewidth',1.5,'color',lcolor(2,:));
        hl = plot(res_box(i).t,res_box(i).Tlower,'linewidth',1.5,'color',lcolor(3,:));
    
        % plot(res_box(i).t,res_box(i).Tupper_min,'linewidth',1.5,'color',lcolor(1,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tupper_max,'linewidth',1.5,'color',lcolor(1,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tinter_min,'linewidth',1.5,'color',lcolor(2,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tinter_max,'linewidth',1.5,'color',lcolor(2,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tlower_min,'linewidth',1.5,'color',lcolor(3,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tlower_max,'linewidth',1.5,'color',lcolor(3,:),'LineStyle',':');
       
        for i_day=1:length(tgt_days)
            xline(tgt_days(i_day),'linestyle','--','color',lcolor(3+i_day,:),'linewidth',2)
        end
    end
    ylim([-2 5.5])
    if mod(i,2) > 0, ylabel('Temperature (^oC)'); end
    if i>n_fjord_runs-2 
        xlabel('Model time (days)'); 
    else
        set(gca,'xticklabels',[]);
   end
    set(gca,'fontsize',fsize)
    if i==1 
        string_legend = {"Shelf","Fjord"};
        for i_day=1:length(tgt_days)
            string_legend{end+1} = sprintf("Model_{%s}",name_days{i_day});
        end
        hl1 = legend([hs, hf, hb],string_legend,'fontsize',fsize,'Location','Southeast'); 
        % title(hl1,sprintf('Profiles at day %d\n(10-day avg.)',tgt_day))
        hl1.NumColumns=2;
    end
    if i==1%n_fjords
        hl2 = legend([hu,hi,hl],{"0-50 m","50-250 m","250-500 m"},'fontsize',fsize,'Location','Northeast'); 
        % title(hl2,'Time series')
        hl2.NumColumns=1;
    end
end


% figure(hf_profiles); exportgraphics(gcf,[figs_path,'temp_profiles_',num2str(2015+i_yr),'_n',num2str(n_C0*n_P0*n_M0*n_gamma),'.png'],'Resolution',300)
% figure(hf_series); exportgraphics(gcf,[figs_path,'temp_series_',num2str(2015+i_yr),'_n',num2str(n_C0*n_P0*n_M0*n_gamma),'.png'],'Resolution',300)
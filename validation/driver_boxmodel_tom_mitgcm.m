%% Driver file for simulating the fjords simulated by Tom using MITgcm


%% Initialise all needed paths and variables
run setup_paths

i_yr = 4; % we want only 2019 data for now, because we only have geometry data for glaciers with data for 2019
time_start = datetime(2019,01,01);
time_end   = datetime(2020,02,15);
input_dt   = 30;
model_dt   = 0.1;
tgt_day = 300; % which day of the 400-day run we want to compare

fjord_letters     = {'B','F','G','M'};
fjords_files_list = {'NUU_80','KAK_10','NAA_20','SES_10'};

%% Compiling all "external" data used (Greenland fjord compilation and Cowton et al., 2023)
run compile_process_fjords % requires data-compilation repository

run load_cowton2023_data

%% Selecting which fjords from which compilation we want, and the time period for discharge data

% These are the IDs of the corresponding fjords above in the "fjords_processed" data structure
% useful function to find them:
% plot_fjords_summary(datasets,fjords_map,fjords_compilation); %plt_handles.cb1.Visible = 'off'; plt_handles.cb2.Visible = 'off'; plt_handles.cb3.Visible = 'off';
fjord_ids = [9,23,24,36];

% will look for matching letters in the table to find which fjords we want
ids_cowton_fjords = zeros(size(fjord_letters));
for i=1:height(meta_table)
    for j=1:length(fjord_letters)
        if strcmp(meta_table.Var1(i),fjord_letters{j}), ids_cowton_fjords(j)=i;end
    end
end

tgt_period = 150:250; %time_glaciers > time_start & time_glaciers < time_end;
qsg_all = qsg_glaciers(tgt_period,:);
taxis_qsg = time_glaciers(tgt_period);
clear qsg_glaciers time_glaciers tgt_period % bit of tidying up

qsg_all    = qsg_all(:,ids_cowton_fjords);
meta_table = meta_table(ids_cowton_fjords,:);


%% Reading the MITgcm outputs

if exist('mitgcm',"var"),       clear mitgcm; end
mitgcm(length(fjords_files_list)) = struct('x',[],'y',[],'z',[],'t',[],'Tprofile',[],'Tupper',[],'Tinter',[],'Tlower',[],'Iprofile',[],'Ivolume',[]);
for i=1:length(fjords_files_list)
    mitgcm_out = load([inputs_path,'/MITgcm/',fjords_files_list{i}]);
    icebergs_out = load([inputs_path,'/MITgcm/icebergArea_',fjords_files_list{i}]);
    dx=400;
    dy=250;
    dz=10;

    x_mitgcm=1:dx:size(mitgcm_out.T,1)*dx;
    y_mitgcm=1:dy:size(mitgcm_out.T,2)*dy;
    z_mitgcm=1:dz:size(mitgcm_out.T,3)*dz;
    t_mitgcm=0:10:400;
    j_fjord = 1:144;

    % Using zero as _fillValue is not a good practice, but should not be an issue
    mitgcm_out.T(mitgcm_out.S == 0) = NaN;
    temp_out = squeeze(mean(mitgcm_out.T(j_fjord,:,:,:),[1,2],'omitnan'));

    mitgcm(i).x = x_mitgcm;
    mitgcm(i).y = y_mitgcm;
    mitgcm(i).z = z_mitgcm;
    mitgcm(i).t = t_mitgcm;
    mitgcm(i).Tprofile = temp_out(:,tgt_day/10);
    mitgcm(i).Tupper = mean(temp_out(1:5,:),1,'omitnan');
    mitgcm(i).Tinter = mean(temp_out(5:25,:),1,'omitnan');
    mitgcm(i).Tlower = mean(temp_out(25:50,:),1,'omitnan');
    mitgcm(i).Iprofile = squeeze(mean(icebergs_out.icebergArea,[1,2],'omitnan'))/(dx*dy);
    mitgcm(i).a = trapz(mitgcm(i).z,mitgcm(i).Iprofile);
    % tetrahedral conversion of area to volume, assuming the iceberg concentration integral as the characteristic length scale
    % mitgcm(i).Ivolume = mitgcm(i).a/(6.*sqrt(6)).*icebergs_out.icebergArea; 
    mitgcm(i).Ivolume = dz.*icebergs_out.icebergArea; 

    clear i_fjord temp_out icebergs_out t_mitgcm z_mitgcm y_mitgcm x_mitgcm dx dy dz% bit of tidying up
end
disp('MITgcm data read.')
%% Preparing the model runs
n_fjords = size(meta_table,1);

if exist('fjord_model',"var"),       clear fjord_model; end
fjord_model(n_fjords) = struct("p",[],"a",[],"f",[],"t",[],"m",[],"c",[]);
for i_fjord=1:n_fjords

    i_fjord_ctd = ids_cowton_fjords(i_fjord);
    % time axis for the run
    t = 0:model_dt:400; % 400 days (as in the MITgcm simulations)

    % compile metadata and default parameters
    m.ID   = meta_table.Var1(i_fjord);
    m.name = meta_table.GreenlandicName(i_fjord);
    m.lon  = meta_table.Longitude__E_(i_fjord);
    m.lat  = meta_table.Latitude__N_(i_fjord);
    [p,~]  = get_model_default_parameters();
    
    % obtain CTD data
    ctd.temp  = ctd_data.Tdata{i_yr,i_fjord_ctd};
    ctd.salt  = ctd_data.Sdata{i_yr,i_fjord_ctd};
    ctd.depth = ctd_data.zdata{i_yr,i_fjord_ctd};
    

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
    f.Qsg = [zeros([1,200/model_dt]), qsg_dtmodel, zeros([1,100/model_dt])];
    clear taxis_qsg_julian taxis_qsg_num taxis_qsg_num_dtmodel qsg_dtmodel % bit of tidying up
    % no meltwater
    % f.Qsg=zeros(size(t));
    % p.P0=0;

    % no icebergs
    % p.M0=0; 
    % f.zi = -800:10:0;
    % % f.zi = f.zi';
    % f.D  = zeros(1,length(t));
    % f.xi = (p.nu0/abs(p.zgl))*exp(p.nu0*f.zi'/abs(p.zgl))/(1-exp(-p.nu0));
    % a.I0 = 0*f.zi;

    % obs-based
    % runtime_axis = datasets.opts.time_start:model_dt:datasets.opts.time_end;
    % [f.zi,ice_d,f.xi,a.I0] = get_total_iceberg_discharge(datasets,fjords_compilation(fjord_ids(i_fjord)).glaciers,runtime_axis',p,a);
    % f.D = ice_d(1:length(t)); % we need to interp/extrapolate the 1-year discharge for 400 days
    
    % MITgcm-based (hybrid: using f.xi and f.zi from here and f.D from obs)
    p.icestatic=1; % icebergs do not evolve in MITgcm
    f.xi = flip(mitgcm(i_fjord).Iprofile./mitgcm(i_fjord).a); % get avg. iceberg concentration in each layer and make its integral equal 1
    f.zi = -flip(mitgcm(i_fjord).z);
    f.D = zeros(size(t)); % only consider the initial iceberg concentration in the fjord, and no discharge afterwards
    D0 = sum(mitgcm(i_fjord).Ivolume(:),'omitnan');
    a.I0 = (D0/(p.W*p.L)).*f.xi;

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

    % % Sanity-check plot to make sure the boxes, T,S profiles, Qsg and D all make sense
    % figure('Position',[100 100 800 500]);
    % ints=[0,-cumsum(a.H0)];
    % subplot(1,3,1); hold on; box on;
    % plot(f.Ts(:,1),f.zs);
    % scatter(a.T0,(ints(1:end-1)+ints(2:end))/2,'filled');
    % for j=1:size(ints,2)
    %     yline(ints(j));
    % end
    % yline(p.zgl,'--k','linewidth',0.5)
    % yline(p.silldepth,'--k','linewidth',0.5)
    % ylim([-p.H 0])
    % subplot(1,3,2); hold on; box on;
    % plot(f.Ss(:,1),f.zs);
    % scatter(a.S0,(ints(1:end-1)+ints(2:end))/2,'filled');
    % for j=1:size(ints,2)
    %     yline(ints(j));
    % end
    % yline(p.zgl,'--k','linewidth',0.5)
    % yline(p.silldepth,'--k','linewidth',0.5)
    % ylim([-p.H 0])
    % subplot(1,3,3); hold on; box on; grid on;
    % plot(t,f.Qsg);
    % yyaxis right
    % yline(D0);
    % legend('Subglacial discharge','Solid-ice discharge');
    % %
    % exportgraphics(gcf,[figs_path,'forcing_summary_MITice_fjord',fjord_letters{i_fjord},'.png'],'Resolution',300)IP
    % 

    clear m p a f % bit of tidying up
end
disp('Model inputs processing done.')
%% Run the model

n_fjord_runs = length(fjord_model);
range_C0 = [1e3,5e3,1e4,5e4,1e5];
range_P0 = 5:5:30;
range_M0 = [5e-9,1e-8,2e-8,5e-8];
range_gamma = 0.1:0.05:0.9;
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
        % ensemble(i_fjord,i_C0,i_P0).s.Tfinal = mean(cur_fjord.s.T(:,end-365:end),2); % last 30 days(?)
        % ensemble(i_fjord,i_C0,i_P0).s.Sfinal = mean(cur_fjord.s.S(:,end-365:end),2);
        % ensemble(i_fjord,i_C0,i_P0).s.Hfinal = mean(cur_fjord.s.H(:,end-365:end),2);
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tfinal = mean(cur_fjord.s.T(:,(tgt_day-5:tgt_day+5)/model_dt),2); % 10-day avg centered at the target day (since MITgcm outputs a 10-day average)
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Sfinal = mean(cur_fjord.s.S(:,(tgt_day-5:tgt_day+5)/model_dt),2);
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Hfinal = mean(cur_fjord.s.H(:,(tgt_day-5:tgt_day+5)/model_dt),2);

        Tregular=NaN([length(mitgcm(i_fjord).z),length(cur_fjord.s.t)]);
        for k=1:length(cur_fjord.s.t)
            ints=[0; cumsum(cur_fjord.s.H(:,k))];
            z_box = (ints(1:end-1)+ints(2:end))/2;
            Tregular(:,k) = interp1(z_box,cur_fjord.s.T(:,k),mitgcm(i_fjord).z,'nearest','extrap');
        end
        
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tupper = mean(Tregular(1:5,:),1,'omitnan');   % avg temp of the upper 50 m
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tinter = mean(Tregular(5:25,:),1,'omitnan');  % avg temp of the 50 - 250 m layer
        ensemble(i_fjord,i_C0,i_P0,i_M0,i_gamma).s.Tlower = mean(Tregular(25:50,:),1,'omitnan'); % avg temp of the 250 - 500 m layer
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
save([outs_path,'boxmodel_runs_MITgcm_MITice_comp_n',num2str(n_combinations),'_day',num2str(tgt_day),'_icebergs_cube'],'-v7.3','ensemble','fjord_model');
disp('Outputs saved.')
% load([outs_path,'boxmodel_runs_MITgcm_obsice_comp_n',num2str(n_combinations),'_day',num2str(tgt_day)]);
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
    tf_box_comp = NaN([length(zf_obs),n_fjord_runs,n_C0,n_P0,n_M0,n_gamma]);
    sf_box_comp = NaN(size(tf_box_comp));
    tupper_box_comp = NaN([length(fjord_model(i).t)-1,n_fjord_runs,n_C0,n_P0,n_M0,n_gamma]);
    tinter_box_comp = NaN(size(tupper_box_comp));
    tlower_box_comp = NaN(size(tupper_box_comp));
    for j=1:size(ensemble,2)
    for k=1:size(ensemble,3)
    for l=1:size(ensemble,4)
    for n=1:size(ensemble,5) % we skip 'm' because it is being used for the metadata structure elsewhere in the code
        if ~isempty(ensemble(i,j,k,l,n).s) & ~isnan(ensemble(i,j,k,l,n).s.Hfinal)
            ints=[0; cumsum(ensemble(i,j,k,l,n).s.Hfinal)];
            z_box = (ints(1:end-1)+ints(2:end))/2;
            tf_box = ensemble(i,j,k,l,n).s.Tfinal; 
            sf_box = ensemble(i,j,k,l,n).s.Sfinal; 
            
            tf_box_comp(:,i,j,k,l,n) = interp1(z_box,tf_box,zf_obs,'nearest','extrap');
            sf_box_comp(:,i,j,k,l,n) = interp1(z_box,sf_box,zf_obs,'nearest','extrap');

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

    res_box(i).t  = t(2:end); %ensemble(i,1,1,1,1).t(2:end);
    res_box(i).zf = z_box;
    res_box(i).tf = median(tf_box_comp(:,i,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).tfmin = min(tf_box_comp(:,i,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).tfmax = max(tf_box_comp(:,i,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).sf = median(sf_box_comp(:,i,:,:,:),[3,4,5,6],'omitnan');


    res_box(i).Tupper = median(tupper_box_comp(:,i,:,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).Tupper_min = min(tupper_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).Tupper_max = max(tupper_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');

    res_box(i).Tinter = median(tinter_box_comp(:,i,:,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).Tinter_min = min(tinter_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).Tinter_max = max(tinter_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');

    res_box(i).Tlower = median(tlower_box_comp(:,i,:,:,:,:),[3,4,5,6],'omitnan');
    res_box(i).Tlower_min = min(tlower_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');
    res_box(i).Tlower_max = max(tlower_box_comp(:,i,:,:,:,:),[],[3,4,5,6],'omitnan');

    res_box(i).id = fjord_model(i).m.ID{1};
    res_box(i).name = fjord_model(i).m.name{1};
    res_box(i).n = n_completed;
end
disp('Postprocessing done.')
%% Plotting results
lcolor = lines(3);
fsize=12;
figure('Name','Temperature comparison','Position',[40 40 1200 400*length(fjord_model)]);
tiledlayout(length(fjord_model),2);
for i=1:n_fjord_runs
    
    nexttile; hold on; box on; grid on
    text(0.02,1.05,sprintf("(%s) %s (%.0f km long)",res_box(i).id,res_box(i).name, fjord_model(i).p.L/1e3),'units','normalized','fontsize',fsize)
    text(0.02,0.05,sprintf("n=%d",res_box(i).n),'Units','normalized','FontSize',fsize)

    % figure; hold on
    % y2 = [res_obs(i).zf; flip(res_obs(i).zf)]';
    % inBetween = [res_box(i).tfmin, flip(res_box(i).tfmax)];
    % hp = fill(flip(inBetween,1), flip(-y2,2), lcolor(3,:),'edgecolor','none','facealpha',0.2);

    hs = plot(res_obs(i).ts,-res_obs(i).zs,'linewidth',1.5,'color',lcolor(1,:));
    hf = plot(res_obs(i).tf,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(2,:));
    hb = plot(res_box(i).tf,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:));
    hm = plot(mitgcm(i).Tprofile,-mitgcm(i).z,'linewidth',1.5,'color','k');
    plot(res_box(i).tfmin,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
    plot(res_box(i).tfmax,-res_obs(i).zf,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
    scatter(0,fjord_model(i).p.zgl,40,'v','filled','MarkerFaceColor','black')
    plot([0 0],[-fjord_model(i).p.H fjord_model(i).p.silldepth],'-k','linewidth',2)
    xlabel('Temperature (^oC)'); ylabel('Depth (m)');
    set(gca,'fontsize',fsize)
    xlim([-2 6])
    ylim([-fjord_model(i).p.H 0])

    nexttile; hold on; box on; grid on
    if length(res_box(i).t) == length(res_box(i).Tupper)
        hbl = plot(res_box(i).t,res_box(i).Tupper,'linewidth',1.5,'color','k'); % the first two lines are dummy for the legend entries
        hml = plot(res_box(i).t,res_box(i).Tupper,'linewidth',1.5,'color','k','LineStyle','--');
    
        hu = plot(res_box(i).t,res_box(i).Tupper,'linewidth',1.5,'color',lcolor(1,:));
        hi = plot(res_box(i).t,res_box(i).Tinter,'linewidth',1.5,'color',lcolor(2,:));
        hl = plot(res_box(i).t,res_box(i).Tlower,'linewidth',1.5,'color',lcolor(3,:));
    
        plot(mitgcm(i).t,mitgcm(i).Tupper,'linewidth',1.5,'color',lcolor(1,:),'LineStyle','--');
        plot(mitgcm(i).t,mitgcm(i).Tinter,'linewidth',1.5,'color',lcolor(2,:),'LineStyle','--');
        plot(mitgcm(i).t,mitgcm(i).Tlower,'linewidth',1.5,'color',lcolor(3,:),'LineStyle','--');
        % plot(res_box(i).t,res_box(i).Tupper_min,'linewidth',1.5,'color',lcolor(1,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tupper_max,'linewidth',1.5,'color',lcolor(1,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tinter_min,'linewidth',1.5,'color',lcolor(2,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tinter_max,'linewidth',1.5,'color',lcolor(2,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tlower_min,'linewidth',1.5,'color',lcolor(3,:),'LineStyle',':');
        % plot(res_box(i).t,res_box(i).Tlower_max,'linewidth',1.5,'color',lcolor(3,:),'LineStyle',':');
       
    end
    ylim([-2 5.5])
    ylabel('Temperature (^oC)'); xlabel('Model time (days)');
    set(gca,'fontsize',fsize)
    if i==1 
        hl1 = legend([hs, hf, hb, hm], {"Shelf","Fjord","Box model","MITgcm"},'fontsize',fsize,'Location','Southeast'); 
        title(hl1,sprintf('Profiles at day %d\n(10-day avg.)',tgt_day))
    end
    if i==1%n_fjords
        hl2 = legend([hbl,hml,hu,hi,hl],{"Box model","MITgcm","0-50 m","50-250 m","250-500 m"},'fontsize',fsize,'Location','Northeast'); 
        title(hl2,'Time series')
        hl2.NumColumns=1;
    end
end


% exportgraphics(gcf,[figs_path,'temp_profiles_MITgcm_MITice_comp_n',num2str(n_combinations),'_day',num2str(tgt_day),'_icebergs_cube.png'],'Resolution',300)
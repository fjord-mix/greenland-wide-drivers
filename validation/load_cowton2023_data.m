%% modelled subglacial discharge
qsg_glaciers = load([inputs_path,'/runoff/runoff.txt']);
time_glaciers = datetime(load([inputs_path,'runoff/runoff_tvec.txt']));

%% CTD profiles
ctd_data   = load([inputs_path,'/fjords/NW_fjord_TSz']);
meta_table = readtable([inputs_path,'/fjords/Fjords_table.xlsx'],'Range','A3:I17'); % I had to unmerge some sells from the original table, otherwise Matlab would not get the header names (why were they merged in the first place?!)

%% MITgcm outputs
% simulations for a subset of fjords from Cowton et al., 2023
fjord_ids_mitgcm  = {'B','F','G','M'}; % [2,6,7,13]
fjords_files_list = {'NUU_80','KAK_10','NAA_20','SES_10'};
tgt_day_mitgcm    = 300; % equivalent post-melt peak in the MITgcm experiments

if exist('mitgcm',"var"),       clear mitgcm; end
mitgcm(length(fjords_files_list)) = struct('x',[],'y',[],'z',[],'t',[],'Tprofile',[],'Tupper',[],'Tinter',[],'Tlower',[],'Iprofile',[],'Ivolume',[]);
for i=1:length(fjords_files_list)
    mitgcm_out = load([inputs_path,'/MITgcm/',fjords_files_list{i}]);
    % icebergs_out = load([inputs_path,'/MITgcm/icebergArea_',fjords_files_list{i}]);
    dx=400;
    dy=250;
    dz=10;

    x_mitgcm=1:dx:size(mitgcm_out.T,1)*dx;
    y_mitgcm=1:dy:size(mitgcm_out.T,2)*dy;
    z_mitgcm=1:dz:size(mitgcm_out.T,3)*dz;
    t_mitgcm=0:10:400;

    % Using zero as _fillValue is not a good practice, but should not be an issue
    mitgcm_out.T(mitgcm_out.S == 0) = NaN;
    mitgcm_out.S(mitgcm_out.S == 0) = NaN;

    % find the ice front, which is the last NaN value in the centre line
    land_ocn_transect = isnan(mitgcm_out.S(:,ceil(size(mitgcm_out.S,2)/2),1,1));
    j_cf = find(land_ocn_transect==1,1,'last');
    
    cf_dist      = 16e3;                  % cf_dist is 16km because it's the furthest away any profile is from the ice front
    j_downstream = j_cf+ceil(cf_dist/dy); % j_downstream is "cf_dist" km downstream from the ice front
    j_fjord      = j_cf+1:j_downstream;   % j_fjord is then j_cf:j_downstream

    temp_out = squeeze(mean(mitgcm_out.T(j_fjord,:,:,:),[1,2],'omitnan'));
    salt_out = squeeze(mean(mitgcm_out.S(j_fjord,:,:,:),[1,2],'omitnan'));

    mitgcm(i).id = fjord_ids_mitgcm{i};
    mitgcm(i).x = x_mitgcm;
    mitgcm(i).y = y_mitgcm;
    mitgcm(i).z = z_mitgcm;
    mitgcm(i).t = t_mitgcm;
    mitgcm(i).Tprofile = temp_out(:,tgt_day_mitgcm/10);
    mitgcm(i).Sprofile = salt_out(:,tgt_day_mitgcm/10);
    mitgcm(i).Tseries = temp_out;
    mitgcm(i).Sseries = salt_out;
    mitgcm(i).Tupper = mean(temp_out(1:5,:),1,'omitnan');
    mitgcm(i).Tinter = mean(temp_out(5:25,:),1,'omitnan');
    mitgcm(i).Tlower = mean(temp_out(25:50,:),1,'omitnan');
    mitgcm(i).Supper = mean(salt_out(1:5,:),1,'omitnan');
    mitgcm(i).Sinter = mean(salt_out(5:25,:),1,'omitnan');
    mitgcm(i).Slower = mean(salt_out(25:50,:),1,'omitnan');

    clear i_fjord j_fjord cf_dist j_cf j_downstream temp_out icebergs_out t_mitgcm z_mitgcm y_mitgcm x_mitgcm dx dy dz% bit of tidying up
end
disp('MITgcm data read.')
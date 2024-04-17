% modelled subglacial discharge
qsg_glaciers = load([inputs_path,'/runoff/runoff.txt']);
time_glaciers = datetime(load([inputs_path,'runoff/runoff_tvec.txt']));

% CTD profiles
ctd_data   = load([inputs_path,'/fjords/NW_fjord_TSz']);
meta_table = readtable([inputs_path,'/fjords/Fjords_table.xlsx'],'Range','A3:I17'); % I had to unmerge some sells from the original table, otherwise Matlab would not get the header names (why were they merged in the first place?!)
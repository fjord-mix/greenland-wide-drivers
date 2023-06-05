function plot_fluxes(fjord_model1,fjord_model2)



model_runtime1 = fjord_model1.s.t(1:size(fjord_model1.s.H,2));
runtime_axis = fjord_model1.m.time_axis;
t0 = convertTo(runtime_axis(1),'datenum');
taxis1 = NaT([size(fjord_model1.s.H,2),1]);
for i_time=1:length(taxis1)
    taxis1(i_time) = datetime(t0+model_runtime1(i_time),'ConvertFrom','datenum');
end
n_steps = length(taxis1);

n_layers=size(fjord_model1.s.H,1);
unit_factor=1e-4;

m=n_layers; n=3;
figure('Position',[20 20 900 900],'Name',fjord_model1.m.name); 
for i_layer=1:n_layers
    subplot(m,n,1+(i_layer-1)*3); hold on; box on
    if i_layer==1, title(sprintf('Volume flux (%0.1e*m^3/s)',unit_factor)); end
    plot(taxis1,real(fjord_model1.s.QVg(i_layer,1:n_steps)).*unit_factor,'linewidth',1.5); 
    plot(taxis1,real(fjord_model1.s.QVs(i_layer,1:n_steps)).*unit_factor,'linewidth',1.5); 
    plot(taxis1,real(fjord_model1.s.QVk(i_layer,1:n_steps)).*unit_factor,'linewidth',1.5); 
    plot(taxis1,real(fjord_model1.s.QVb(i_layer,1:n_steps)).*unit_factor,'linewidth',1.5); 
    line(taxis1,0*model_runtime1,'linestyle',':');
    % icebergs dont have a volume flux
    text(0.02,0.9,sprintf('Layer %d',i_layer),'Units','normalized')    
end


for i_layer=1:n_layers
    subplot(m,n,2+(i_layer-1)*3); hold on; box on
    if i_layer==1, title(sprintf('Temperature flux (%0.1e*???)',unit_factor)); end
    plot(taxis1,real(fjord_model1.s.QTg(i_layer,1:n_steps)).*unit_factor,'linewidth',1.5); 
    plot(taxis1,real(fjord_model1.s.QTs(i_layer,1:n_steps)).*unit_factor,'linewidth',1.5); 
    plot(taxis1,real(fjord_model1.s.QTk(i_layer,1:n_steps)).*unit_factor,'linewidth',1.5); 
    plot(taxis1,real(fjord_model1.s.QTb(i_layer,1:n_steps)).*unit_factor,'linewidth',1.5);     
    plot(taxis1,real(fjord_model1.s.QTi(i_layer,1:n_steps)).*unit_factor,'linewidth',1.5); 
    line(taxis1,0*model_runtime1,'linestyle',':');
    text(0.02,0.9,sprintf('Layer %d',i_layer),'Units','normalized')
end


for i_layer=1:n_layers
    subplot(m,n,3+(i_layer-1)*3); hold on; box on
    if i_layer==1, title(sprintf('Salinity flux (%0.1e*???)',unit_factor)); end
    plot(taxis1,real(fjord_model1.s.QSg(i_layer,1:n_steps)).*unit_factor,'linewidth',1.5); 
    plot(taxis1,real(fjord_model1.s.QSs(i_layer,1:n_steps)).*unit_factor,'linewidth',1.5); 
    plot(taxis1,real(fjord_model1.s.QSk(i_layer,1:n_steps)).*unit_factor,'linewidth',1.5); 
    plot(taxis1,real(fjord_model1.s.QSb(i_layer,1:n_steps)).*unit_factor,'linewidth',1.5);     
    plot(taxis1,real(fjord_model1.s.QSi(i_layer,1:n_steps)).*unit_factor,'linewidth',1.5); 
    line(taxis1,0*model_runtime1,'linestyle',':');
    text(0.02,0.9,sprintf('Layer %d',i_layer),'Units','normalized')
end

legend('plume','shelf','mixing','artificial','icebergs','Location','east');



% subplot(m,n,2); plot(taxis1,fjord_model1.s.QTb,'linewidth',1.5); hold on;
% ylabel('Temperature flux (??)'); 
% subplot(m,n,3); plot(taxis1,fjord_model1.s.QSb,'linewidth',1.5); hold on;
% ylabel('Salinity flux (??)'); xlabel('Model time (days)')


end
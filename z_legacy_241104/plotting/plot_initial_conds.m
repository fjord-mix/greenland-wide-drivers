function hf = plot_initial_conds(hf,fjord_model)

t0 = fjord_model.m.time_axis(1);
T0_box=flip(fjord_model.a.T0);
S0_box=flip(fjord_model.a.S0);
H0_box=fjord_model.a.H0;
Z_box = flip(-cumsum(H0_box));

T_sh_profile = fjord_model.f.Ts(:,1);
S_sh_profile = fjord_model.f.Ss(:,1);
Z_sh_profile = fjord_model.f.zs;

% Copying box values to their respective depths to compare with shelf profile
% NOTE: we cannot really interpolate here, otherwise we won't see the REAL values!
T0_profile=Z_sh_profile.*0;
S0_profile=Z_sh_profile.*0;
for i=1:length(Z_box)
    T0_profile(Z_sh_profile<Z_box(i) & T0_profile==0) = T0_box(i);
    S0_profile(Z_sh_profile<Z_box(i) & S0_profile==0) = S0_box(i);
end
T0_profile(T0_profile==0) = T0_box(end);
S0_profile(S0_profile==0) = S0_box(end);

depths_not_in_fjord=Z_sh_profile<Z_box(1);
T_sh_profile(depths_not_in_fjord)=[];
S_sh_profile(depths_not_in_fjord)=[];
Z_sh_profile(depths_not_in_fjord)=[];
T0_profile(depths_not_in_fjord)=[];
S0_profile(depths_not_in_fjord)=[];


if isempty(hf)
    hf = figure('Name',fjord_model.m.name);
else
    set(gcf,hf)
end

subplot(1,2,1); hold on; box on;
% yline(Z_box,':k','linewidth',0.5);
hs=plot(T_sh_profile,Z_sh_profile,'-b','linewidth',1.5);
hm=plot(T0_profile,Z_sh_profile,':b','linewidth',1.5);
% text(0.95,0.95,date2str(t0),'Units','normalized')
ylabel('Depth (m)');
xlabel('Temperature (^oC)')

subplot(1,2,2); hold on; box on;
% yline(Z_box,':k','linewidth',0.5);
plot(S_sh_profile,Z_sh_profile,'-r','linewidth',1.5);
plot(S0_profile,Z_sh_profile,':r','linewidth',1.5);
ylabel('Depth (m)');
xlabel('Salinity (-)')
legend([hs,hm],{'shelf','model'},'location','southwest');

end
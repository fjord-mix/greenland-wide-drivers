function [hf_t,hf_s] = plot_ens_mitgcm_boxed(res_box,fjord_model,n_fjord_runs)

    fsize=14;
    %% Plotting temperature
    hf_t = figure('Position',[40 40 700 600]);
    ht_t = tiledlayout('flow');
    hb = plot(res_box(1).t,res_box(1).Tbox(1,:),'-k');
    hm = plot(res_box(1).t,res_box(1).Tmitgcm(1,:),'--k');
    for i_fjord=1:n_fjord_runs
        nexttile; hold on; box on;
        text(0.02,1.05,sprintf("(%s) %s (%.0f km long)",res_box(i_fjord).id,res_box(i_fjord).name, fjord_model(i_fjord).p.L/1e3),'units','normalized','fontsize',fsize)
        
        plot(res_box(i_fjord).t,res_box(i_fjord).Tbox,'linewidth',1.5);
        set(gca,'ColorOrderIndex',1)
        plot(res_box(i_fjord).t,res_box(i_fjord).Tmitgcm,'linewidth',1.5,'linestyle','--')
        set(gca,'fontsize',fsize)
        if i_fjord==1
            legend('Layer 1','Layer 2','Layer 3','fontsize',fsize)
        end
        xlabel('Model time (days)')
        ylabel('Temperature (^oC)')
    end
    % handle_models = [hb; hm];
    % lbl_models = {'Box model','MITgcm'};
    % legend(handle_models,{'Box model','MITgcm'})
    % ht.TileSpacing='compact';
    % ht.Padding='compact';

    %% Plotting salinity
    hf_s = figure('Position',[40 40 700 600]);
    ht_s = tiledlayout('flow');
    hb = plot(res_box(1).t,res_box(1).Tbox(1,:),'-k');
    hm = plot(res_box(1).t,res_box(1).Tmitgcm(1,:),'--k');
    for i_fjord=1:n_fjord_runs
        nexttile; hold on; box on;
        text(0.02,1.05,sprintf("(%s) %s (%.0f km long)",res_box(i_fjord).id,res_box(i_fjord).name, fjord_model(i_fjord).p.L/1e3),'units','normalized','fontsize',fsize)
        
        plot(res_box(i_fjord).t,res_box(i_fjord).Sbox,'linewidth',1.5);
        set(gca,'ColorOrderIndex',1)
        plot(res_box(i_fjord).t,res_box(i_fjord).Smitgcm,'linewidth',1.5,'linestyle','--')
        set(gca,'fontsize',fsize)
        if i_fjord==1
            legend('Layer 1','Layer 2','Layer 3','fontsize',fsize)
        end
        xlabel('Model time (days)')
        ylabel('Salinity')
    end

end
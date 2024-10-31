function [P, weightY, weightX] = CCA_plots(Project, mode2plot, weighttype, savepath, beh_labels)

    % [P, weightY, weightX] = CCA_plots('_Project_BodyWandering', 1, 'correlation', '/home/leah/Git/BodyWanderingCCA/figures/CCA_BWitems/', {'Future', 'Past', 'Self', 'Other', 'Pos', 'Neg', 'Words', 'Vivid', 'Vague', 'Focus', 'Ruminate', 'Distant', 'Image', 'Arousal', 'Body', 'Breath', 'Heart', 'Movement', 'Bladder', 'Skin', 'Stomach'}) 

    %clear all
    Project = '_Project_BodyWandering';
    mode2plot = 1; 
    weighttype = 'correlation';  % 'weight'; 
    savepath = '/home/leah/Git/BodyWanderingCCA/figures/CCA_BWitems/'
    beh_labels = {'Future', 'Past', 'Self', 'Other', 'Pos', 'Neg', 'Words', 'Vivid', 'Vague', 'Focus', 'Ruminate', 'Distant', 'Image', 'Arousal', 'Body', 'Breath', 'Heart', 'Movement', 'Bladder', 'Skin', 'Stomach'}; 
   
    %% plotting
    % Set path for plotting and the BrainNet Viewer toolbox
    set_path('plot', 'brainnet');
    addpath('/Users/au704655/Documents/Packages/spm12/')
    
    % % load cfg file for effect of interest
    filepath = ['/home/leah/Git/connectivity_tools/toolboxes/matlab/cca_pls_toolkit-master/' Project '/framework/'];
    framework = 'cca_pca_holdout5-0.40_subsamp5-0.40'; 
    cfg = load([filepath, framework, '/cfg_1.mat']);
    
    % Load res
    res.dir.frwork = cfg.cfg.dir.frwork;
    res.frwork.level = mode2plot;
    % get the correlation between the input variables and the latent variables/projections
    res.gen.weight.type = weighttype; % 'correlation' for structure correlations/loadings & 'weight' for true model weights
    % res.gen.selectfile = 'none';
    % res.gen.weight.flip = 1;
    res = res_defaults(res, 'load');
    
    %% Plot scatterplot of data projections (Connectome * Factor Canonical Variates)
    % P = projections 
    P = plot_proj(res, {'X' 'Y'}, res.frwork.level, 'osplit', res.frwork.split.best, ...
        'training+test', '2d_group', 'gen.axes.FontSize', 20, ...
        'gen.legend.FontSize', 20, 'gen.legend.Location', 'NorthWest', ... 
        'proj.scatter.SizeData', 120, 'proj.scatter.MarkerEdgeColor', 'k', ...
        'proj.scatter.MarkerFaceColor',  [0.3 0.3 0.9; 0.9 0.3 0.3], ... %[[224, 206, 253]/255; [87, 87, 241]/255], ... %'#c8b6ff' , '#957fef'], ... %[0.3 0.3 0.9; 0.9 0.3 0.3] % [242/255 191/255 211/255; 173/255 31/255 87/255
        'proj.xlabel', 'Connectivity Variate', ...  %'Brain Variate'
        'proj.ylabel', 'Body-Wandering Variate', ... %'Beh Variate'
        'proj.lsline', 'on');
        load([res.dir.frwork, sprintf('/res/level%d/model_1.mat', mode2plot)])
        title({sprintf('r = %.03f , p = %.03f', round(correl(res.frwork.split.best),3), round(res.stat.pval(res.frwork.split.best),3)),},'fontsize',18)
        
        xlim([min(P(:,1))-0.02 max(P(:,1))+0.02])
        ylim([min(P(:,2))-0.02 max(P(:,2))+0.02])
    
    saveas(gcf,[savepath, 'scatter_CCAvariates_testtrain.png'])
    
    %    % Plot data projections across levels (and averaged across modalities 
    %    % in a level after standardization)
    %    plot_proj(res, {'X' 'Y'; 'X' 'Y'}, [1 1; 2 2], 'osplit', res.frwork.split.best, ...
    %              'none', '2d', 'proj.multi_label', 1);
    
    %% Plot behavioural weights as vertical bar plot (function from cca-pls toolkit - alternative barplot later)
    %weightY = plot_weight(res, 'Y', 'behav', res.frwork.split.best, 'behav_vert', ...
    %'gen.axes.FontSize', 20, 'gen.legend.FontSize', 12); % ...
    %'behav.weight.sorttype', 'sign', 'behav.weight.numtop', 20);
    %'gen.axes.YLim', [-0.004 0.013], ...
    %if strcmp(weighttype, 'correlation')
    %    title('Behavioural CCA loadings (structure correlations)')
    %    saveas(gcf, [savepath, 'bar_CCAloadings_beh.png'])
    %elseif strcmp(weighttype, 'weight')
    %    title('Behavioural CCA weights')
    %    saveas(gcf, [savepath, 'bar_CCAweights_beh.png'])
    %end

    %% calculate weights
    weightX = calc_weights('X', res.frwork.split.best, res, cfg.cfg);
    weightY = calc_weights('Y', res.frwork.split.best, res, cfg.cfg);
   
    %% Plot scatterplot of Connectome * Mind-Wanding Canonical Variates and color by most influential variable
    % Edit figure properties
    fontName    = 'Helvetica';
    fontSize    = 18;
    dotsize         = 65;
    
    plotX           = P(:,1);%vargout.Q'*vargout.U(:,mode2plot); % connectome projections
    plotY           = P(:,2);%vargout.Q'*vargout.V(:,mode2plot); % body/mind wandering items projections
    rawdata_path = ['/home/leah/Git/connectivity_tools/toolboxes/matlab/cca_pls_toolkit-master/', Project, '/data/'];
    peak_idx = find(abs(weightY) == max(abs(weightY)));
    peak_var      = load([rawdata_path, 'Y.mat']); peak_var = peak_var.Y(:,peak_idx);
    
    figure
    set(gcf,'color','w');
    scatter(plotX, plotY, dotsize, peak_var,'filled') %%%%%%%%%%%%%%%%%%
    
    load([filepath, framework, sprintf('/res/level%d/model_1.mat',mode2plot)])
    title({sprintf('CCA Mode %d', mode2plot), sprintf('r = %.03f , p = %.03f', round(correl(res.frwork.split.best),3), round(res.stat.pval(res.frwork.split.best),3))},'fontsize',fontSize)

    ylabel(sprintf('Body-Wandering Variate'),'fontsize',fontSize); %'Beh Variate' 
    xlabel(sprintf('Connectivity Variate'),'fontsize',fontSize); %'Brain Variate' 
    set(gca,'FontName',fontName,'fontsize',fontSize)
    legend(beh_labels{peak_idx},'Location','southwest');legend boxoff ; %%
    
    xlim([min(P(:,1))-0.02 max(P(:,1))+0.02])
    ylim([min(P(:,2))-0.02 max(P(:,2))+0.02])
     
    saveas(gcf, [savepath, sprintf('scatter_CCAvariates_%s.png',beh_labels{peak_idx})])
    
    
    %% Heatmap of behavioural weights
    % Edit figure properties
    fontName    = 'Helvetica';
    fontSize    = 18;
    
    subplot(1,2,1);
    Aload = weightX;%wX_corr;%corr(vargout.Q'*vargout.Y0,vargout.Q'*vargout.U);
    imagesc(Aload); colorbar %imagesc(Aload(res.frwork.split.best,:)'); colorbar
    ylabel(sprintf('Original variables'),'fontsize',fontSize);
    xlabel(sprintf('Canonical variables'),'fontsize',fontSize);
    title({'Imaging side',''}','fontsize',fontSize);
    % daspect([1 1 1])
    
    set(gca,'FontName',fontName,'fontsize',fontSize-8)
    
    subplot(1,2,2)
    Bload = weightY;%wY_corr; %corr(vargout.Q'*vargout.X0,vargout.Q'*vargout.V);
    imagesc(Bload); colorbar %imagesc(Bload(res.frwork.split.best,:)'); colorbar
    title({'Task side',''},'fontsize',fontSize-2);
    ylabel(sprintf('Original variables '),'fontsize',fontSize-2);
    xlabel(sprintf('Canonical variables '),'fontsize',fontSize-2);
    %daspect([1 1 1])

    % Set color axis limits to ensure symmetry
    caxis([-max(abs(Bload(:))), max(abs(Bload(:)))]);

    varNames = beh_labels;

    yticks(1:length(varNames))
    yticklabels(varNames)
    %ytickangle(45)
    set(gcf,'color','w');
    set(gca,'FontName',fontName,'fontsize',fontSize-10)
    
    if strcmp(weighttype, 'correlation')
        sgtitle('CCA Loadings (structure correlations)','fontsize',fontSize+2) 
        saveas(gcf, [savepath, 'heatmap_CCAloadings_equalcolorbar.png'])
        save([savepath, 'brainloadings'], 'weightX')
        save([savepath, 'psychloadings'], 'weightY')
    elseif strcmp(weighttype, 'weight')
        sgtitle(sprintf(sprintf('CCA %ss', weighttype)),'fontsize',fontSize+2) 
        saveas(gcf, [savepath, 'heatmap_CCAweights_equalcolorbar.png'])
        save([savepath, 'brainweights'], 'weightX')
        save([savepath, 'psychweights'], 'weightY')
    end

    % save projections/variate also
    save([savepath, 'CCAvariate'], 'P')


    %% ordered bar plot of beh loadings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~, I] = sort(Bload, 'descend');
    %I = I([1:5, 31:37]); % if want to select subset to plot
    figure;
    bh = bar(reordercats(categorical(varNames(I)), varNames(I)), Bload(I), 'FaceColor', 'flat');

    % MAKE a custom colormap of your choice
    % Define your two RGB color codes
    color1_rgb = [255, 144, 179]/255; % Example: Red
    color2_rgb = [110, 68, 255]/255; % Example: Blue

    % Number of colors in the colormap
    num_colors = length(I); % Assuming 'I' is defined in your context

    % Generate the custom colormap
    custom_colormap = flipud([linspace(color1_rgb(1), color2_rgb(1), num_colors)', ...
                  linspace(color1_rgb(2), color2_rgb(2),  num_colors)', ...
                  linspace(color1_rgb(3), color2_rgb(3), num_colors)']);

    % Use this custom colormap
    colormap(custom_colormap);

    % Define a custom colormap of your choice
    %custom_colormap = flipud(cool(length(I)));  % flipud reverse colourmap %%parula
    %colormap(custom_colormap);


    % Set the colors of each bar individually based on Bload values
    for i = 1:length(I)
        bh.CData(i, :) = custom_colormap(i, :);
    end

    ax = gca; % Get current axes
    %ax.XAxis.FontSize = 14; 
    ax.FontSize = 14; % This will affect both axes and the title


    title('CCA Loadings (structure correlations)');
    ylabel('CCA loadings');
    xlabel('Body/Mind-Wandering loadings'); 
    box off;

    if strcmp(weighttype, 'correlation')
        saveas(gcf, [savepath, 'barplot_CCAloadings_cool.png']);
    elseif strcmp(weighttype, 'weight')
        saveas(gcf, [savepath, 'barplot_CCAweights_cool.png']);
    end


    %% ordered vertical bar plot (coloured bars by intensity) %%%%%%%%%%%%%%%%
    [~, I] = sort(Bload, 'descend');
    figure;
    bh = barh(reordercats(categorical(varNames(I)), varNames(I)), Bload(I), 'FaceColor', 'flat');

    set(gcf, 'Position', [100, 100, 600, 800]);  % Adjust the values as needed
    set(gca, 'Position', [0.7, 0.1, 0.25, 0.8]);

    xlabel('CCA loadings');
    ylabel('Body/Mind-Wandering variables'); 
    box off;

    custom_colormap = flipud(cool(length(I)));  %flipud reverses colours %%flipud(cool(length(I)));
    colormap(custom_colormap);
    % Set the colors of each bar individually based on Bload values
    for i = 1:length(I)
        bh.CData(i, :) = custom_colormap(i, :);
    end

    % Save the figure
    if strcmp(weighttype, 'correlation')
        saveas(gcf, [savepath, 'barplot_CCAloadings_vert_cool.png']);
    elseif strcmp(weighttype, 'weight')
        saveas(gcf, [savepath, 'barplot_CCAweights_vert_cool.png']);
    end

    
    %% Plot hyperparameter surface for grid search results
    plot_paropt(res, 1, {'correl', 'simwx', 'simwy'}, ...
    'gen.figure.Position', [500 600 1200 400], 'gen.axes.FontSize', 20, ...
    'gen.axes.XScale', 'log', 'gen.axes.YScale', 'log');
    saveas(gcf, [savepath, 'hyperparameters.png']);
    
    % plotting grid search results
    %plot_paropt(res, 1, {'trcorrel', 'correl'}, ...
    %'gen.figure.Position', [500 600 800 400], 'gen.axes.FontSize', 20, ...
    %'gen.axes.XScale', 'log', 'gen.axes.YScale', 'log');


    %% wordcloud of beh loadings
    %behloadings = load("/home/leah/Git/BodyWanderingCCA/figures/CCA_BWitems/psychloadings.mat");
    %behlabels = beh_labels;
    %%%%cmap = cool(length(behloadings.weightY));
    %%%%cl = interp1(linspace(min(abs(behloadings.weightY)), max(abs(behloadings.weightY)), length(cmap)), cmap, abs(behloadings.weightY), 'nearest');
    %cl = flipud(custom_colormap)
    %figure;
    %wordcloud(behlabels, abs(behloadings.weightY), 'Color', cl); 
    %saveas(gcf, [savepath, 'behloadings_WordCloud.png']);


    %% Get CCA results for table
    res = res_defaults(res, 'load');

    % always look at most sig data split (res.frwork.split.best)? 
    load([res.dir.frwork, sprintf('/res/level%d/model_1.mat', mode2plot)])
    
    % Weight Stability - Connecitivity:
    result1 = sprintf('Brain (X) Weight Stability = %0.3f', round(mean(simwx(res.frwork.split.best,:)),3))
    
    % Explained Variance - Connectivity (in-sample - recommended):
    result2 = sprintf('Brain (X) Explained Variance (in-sample) = %0.3f', round(trexvarx(res.frwork.split.best),3))
    
    % Weight Stability - Mind-Wandering:
    result3 = sprintf('Behaviour (Y) Weight Stability = %0.3f', round(mean(simwy(res.frwork.split.best,:)),3))
    
    % Explained Variance - Mind-Wandering (in-sample - recommended):
    result4 = sprintf('Behaviour (Y) Explained Variance (in-sample) = %0.3f', round(trexvary(res.frwork.split.best),3))
    
    % In-Sample Correlation
    result5 = sprintf('CCA in-sample Correlation = %0.3f', round(trcorrel(res.frwork.split.best),3))
    % no p-val for in-sample correlation (because permutations based on predictive framework not descriptive)
    
    % Out-of-Sample Correlation - model generalizability 
    result6 = sprintf('CCA out-of-sample Correlation = %0.3f, p = %0.3f', round(correl(res.frwork.split.best),3), round(res.stat.pval(res.frwork.split.best),3))
    
    % Robustness (number of significant data splits) of the CCA model
    result7 = sprintf('Robustness of CCA model = %d percent', (sum(res.stat.pval < (0.05 / length(res.stat.pval))) / length(res.stat.pval)) * 100)

    % Open a text file for writing
    fileID = fopen([savepath, 'CCAmode_results.txt'], 'w');
    % Write the results to the file
    fprintf(fileID, '%s\n%s\n%s\n%s\n%s\n%s\n%s\n', result1, result2, result3, result4, result5, result6, result7);
    fclose(fileID); % Close the file


end

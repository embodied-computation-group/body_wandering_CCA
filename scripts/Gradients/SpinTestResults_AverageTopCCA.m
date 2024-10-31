%% load results from spin test

% load gradient averages and spin-test null distribution results 
savepath = '/Users/au704655/Documents/Body_wandering/Figures/Gradients/Gradient1/';
load('/Users/au704655/Documents/Body_wandering/Data/gradients/gradientcort_surface/corticalsubcortical/gradient1/gradient1_spin.mat');

% load top 1% of parcel regions from CCA
top_parcels = readmatrix('/Users/au704655/Documents/Body_wandering/Data/gradients/gradient_averages_cca.csv'); top_parcels = top_parcels(:,1); % load just for parcel numbers

% select gradient spin-test results for top 1% CCA parcels only (average gradient values and null-distribution averages)
gradients_averages = gradient_spin.average_gradients(ismember(gradient_spin.parcelnums, top_parcels),2);
gradients_nullaverages = gradient_spin.average_gradients_perm(ismember(gradient_spin.parcelnums, top_parcels),:);

% parcel labels 
parcel_labels = readtable('/Users/au704655/Documents/Body_wandering/Data/gradients/gradient_averages_cca.csv'); parcel_labels = parcel_labels(:,[1,7]); % load just for parcel labels
missing_parcels = setdiff(top_parcels, gradient_spin.parcelnums); % caudate L/R missing in surface space(7001/7005) 
parcel_labels = parcel_labels(~ismember(double(parcel_labels{:,1}), missing_parcels), :);
parcel_labels_all = readtable('/Users/au704655/Documents/Body_wandering/Data/gradients/abend216_yeo7+InteroNetwork.csv');
parcel_labels_all = parcel_labels_all(ismember(double(parcel_labels_all{:,5}), double(parcel_labels{:,1})), :);

%% stats
grandaverage = mean(gradients_averages); % average gradient-values across top 1% CCA regions
grandaverage_null = mean(gradients_nullaverages, 1); % same average as above for 2500 null average distribution
p_value = sum(abs(grandaverage_null) >= abs(grandaverage)) / length(grandaverage_null)

% % 5 p-values correct for multiple comparisons
% p_values_all = [0, 0.0480, 0.2416, 0.2052, 0.1888];
% % bonferroni
% p_values_all < (0.05/length(p_values_all))
% % FDR Benjamini-Hochberg
% p_value_all_fdr = mafdr(p_values_all, 'BHFDR', true);
% p_value_all_fdr < 0.05

%% plots %%
% save plots (null distr, bar/box-plots, parcel plots)

%% plot gradient average as compared to spin-test null averages distribution
if p_value < 0.05
    % Extract the null distribution for the chosen gradient
    null_distribution = grandaverage_null;
    % Extract the observed gradient average
    gradient_value = grandaverage;
    % Create a histogram of the null distribution
    figure;
    h = histogram(null_distribution, 30); % 30 bins for the histogram
    h.FaceColor = [0 0 0.5]; 
    hold on;
    y_limits = ylim;
    plot([gradient_value, gradient_value], y_limits, 'r', 'LineWidth', 2);
    % Add labels and title
    xlabel_handle = xlabel('Spin-Test Null Gradient Averages');
    set(xlabel_handle, 'Units', 'normalized', 'Position', [0.5, -0.075, 0]);
    ylabel_handle = ylabel('Frequency');
    set(ylabel_handle, 'Units', 'normalized', 'Position', [-0.09, 0.5, 0]);
    title_handle = title('Gradient-1 Spin-Test Averages of Top CCA Regions');
    set(title_handle, 'Units', 'normalized', 'Position', [0.5, 1.025, 0]);
    %legend('Null Distribution', 'Observed Gradient Value');
    set(gca, 'FontSize', 15);
    hold off;
    box off
    saveas(gcf, sprintf('%sspintest_nulldistr_TopCCAAverage.png', savepath));
end

%% significant parcel surface brain plot
%Set up paths
fshome = '/Applications/freesurfer/7.4.1/'; %getenv('FREESURFER_HOME');
fsmatlab = sprintf('%s/matlab',fshome);
path(path,fsmatlab);
addpath('/Users/au704655/Documents/Body_wandering/Scripts/spin-test-master/scripts')

% gradient atlas in surface space
datal = load_mgh('/Users/au704655/Documents/Body_wandering/Data/gradients/gradientcort_surface/corticalsubcortical/gradient1/gradient1cort_surface.L.mgh'); 
datar =  load_mgh('/Users/au704655/Documents/Body_wandering/Data/gradients/gradientcort_surface/corticalsubcortical/gradient1/gradient1cort_surface.R.mgh');
% schaefer-subcort combined atlas in surface space
lh_parcels = load_mgh('/Users/au704655/Documents/Body_wandering/Data/gradients/schaefer200subcort_atlas.L.mgh');
rh_parcels = load_mgh('/Users/au704655/Documents/Body_wandering/Data/gradients/schaefer200subcort_atlas.R.mgh');

% create mask of top 1% CCA regions
lh_mask = zeros(length(datal),1); %datal; %
rh_mask = zeros(length(datal),1); %datar; %
parcels2plot = table2array(parcel_labels(:,1));
for i = 1:length(parcels2plot)
    if any(parcels2plot(i) == unique(lh_parcels))
        lh_mask(lh_parcels == parcels2plot(i)) = datal(lh_parcels == parcels2plot(i)); 
    elseif any(parcels2plot(i) == unique(rh_parcels))
        rh_mask(rh_parcels == parcels2plot(i)) = datar(rh_parcels == parcels2plot(i)); 
    end
end

%extract the vertice coordinates from the corresponding surface for freesurfer plotting
fshome = '/Applications/freesurfer/7.4.1/'; %getenv('FREESURFER_HOME');
[verticespl, facespl] = freesurfer_read_surf(fullfile(fshome,'subjects/fsaverage/surf/lh.pial'));
[verticespr, facespr] = freesurfer_read_surf(fullfile(fshome,'subjects/fsaverage/surf/rh.pial'));

% custommap=colormap('jet');
% subplot(2,2,1);
% mincol=min(lh_mask);
% maxcol=max(lh_mask);
% plotFSsurf(facespl,verticespl,lh_mask,custommap,mincol,maxcol,[-90 0]);
% title('Significant Gradients Left');
% subplot(2,2,2);
% plotFSsurf(facespl,verticespl,lh_mask,custommap,mincol,maxcol,[90 0]);
% title('Significant Gradients Left');
% 
% subplot(2,2,3);
% mincol=min(rh_mask);
% maxcol=max(rh_mask);
% plotFSsurf(facespl,verticespl,rh_mask,custommap,mincol,maxcol,[-90 0]);
% title('Significant Gradients Right');
% subplot(2,2,4);
% plotFSsurf(facespl,verticespl,rh_mask,custommap,mincol,maxcol,[90 0]);
% title('Significant Gradients Right');

%% gray brain background
% Define the colormap and its size
colormap_name = 'parula'; % You can choose any colormap
cmap = colormap(colormap_name);
num_colors = size(cmap, 1);

% Define subplot positions for the plots
subplot_positions = [...
    0.1 0.55 0.35 0.35;   % Left Hemisphere - view(90, 0)
    0.55 0.55 0.35 0.35;  % Left Hemisphere - view(-90, 0)
    0.1 0.1 0.35 0.35;    % Right Hemisphere - view(90, 0)
    0.55 0.1 0.35 0.35];  % Right Hemisphere - view(-90, 0)

% Data to plot
masks = {lh_mask, lh_mask, rh_mask, rh_mask}; % Cell array of masks
vertices = {verticespl, verticespl, verticespr, verticespr}; % Cell array of vertex data
faces = {facespl, facespl, facespr, facespr}; % Cell array of face data
views = {[90, 0], [-90, 0], [90, 0], [-90, 0]}; % Views for each subplot

% Loop through each subplot
for k = 1:4
    % Create the subplot
    axes('Position', subplot_positions(k, :));
    
    % Get the current mask and vertices
    mask = masks{k};
    vertex_data = vertices{k};
    
    % Normalize mask values to colormap index range [1, num_colors]
    mask_normalized = (mask - min(mask(mask ~= 0))) / ...
                      (max(mask(mask ~= 0)) - min(mask(mask ~= 0))) * ...
                      (num_colors - 1) + 1;

    % Initialize vertex colors with gray
    vertex_colors = repmat([0.8 0.8 0.8], size(vertex_data, 1), 1);

    % Map non-zero mask values to colors using the colormap
    for i = 1:size(vertex_data, 1)
        if mask(i) ~= 0
            index = round(mask_normalized(i)); % Convert normalized value to colormap index
            vertex_colors(i, :) = cmap(index, :); % Set RGB color from colormap
        end
    end

    % Plot the surface
    trisurf(faces{k}, vertex_data(:,1), vertex_data(:,2), vertex_data(:,3), ...
            'FaceVertexCData', vertex_colors, 'FaceColor', 'interp', 'EdgeColor', 'none');
    
    % Adjust colormap and color axis
    colormap(cmap);
    caxis([min(mask(mask ~= 0)) max(mask(mask ~= 0))]);

    % Add a colorbar
    cbar = colorbar;
    set(cbar, 'Location', 'southoutside');

    % Set additional plot settings
    axis equal;
    view(views{k}(1), views{k}(2));
    lighting phong;
    camlight;
    axis off;
    set(gca, 'Box', 'off');
    
    % Set title based on the hemisphere
    if k <= 2
        title('Left Hemisphere');
    else
        title('Right Hemisphere');
    end
end

saveas(gcf, sprintf('%sTopCCAAverage_surfacebrainplot.png', savepath));


%% boxplot of averages of Top 1% CCA regions
gradients = 1:5;
% load top 1% of parcel regions from CCA
top_parcels = readmatrix('/Users/au704655/Documents/Body_wandering/Data/gradients/gradient_averages_cca.csv'); top_parcels = top_parcels(:,1);
for i = 1:length(gradients)
    load(sprintf('/Users/au704655/Documents/Body_wandering/Data/gradients/gradientcort_surface/corticalsubcortical/gradient%d/gradient%d_spin.mat', gradients(i), gradients(i))); 
   % select gradient spin-test results for top 1% CCA parcels only
    gradients_averages_all(:,i) = gradient_spin.average_gradients(ismember(gradient_spin.parcelnums, top_parcels),2);
end
% Create the box plot for each gradient
figure;
h = boxplot(gradients_averages_all, 'Labels', {'Gradient 1', 'Gradient 2', 'Gradient 3', 'Gradient 4', 'Gradient 5'}, ...
    'Notch', 'on', 'Colors', [0 0 0.5], 'Symbol', '');
% Customize the appearance
box off; % Remove the box around the plot
set(gca, 'FontSize', 15); % Increase font size for better readability
% Change box colors and thicken lines
h1 = findobj(gca, 'Tag', 'Box');
for j = 1:length(h1)
    patch(get(h1(j), 'XData'), get(h1(j), 'YData'), [0 0 0.5], 'FaceAlpha', .5);
end
set(h, 'LineWidth', 1.5); % Set the line width of the boxes
% Label the axes and title
ylabel_handle = ylabel('Average Gradient Values');
title_handle = title('Gradient Averages of Top CCA Parcels');
set(ylabel_handle, 'Units', 'normalized', 'Position', [-0.075, 0.5, 0]);
set(title_handle, 'Units', 'normalized', 'Position', [0.5, 1.025, 0]);

% Adjust the position of x-axis labels for spacing
xticks = get(gca, 'XTick'); % Get current x-tick positions
xtick_labels = get(gca, 'XTickLabel'); % Get current x-tick labels
% Clear the default x-tick labels
set(gca, 'XTickLabel', []);
% Calculate the spacing value
y_limits = get(gca, 'YLim'); % Get y-axis limits
spacing = y_limits(1) - 0.025 * (y_limits(2) - y_limits(1)); % Set spacing below the x-axis
% Manually place x-tick labels with increased spacing
for i = 1:length(xticks)
    text(xticks(i), spacing, xtick_labels{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 15);
end

saveas(gcf, '/Users/au704655/Documents/Body_wandering/Figures/Gradients/Gradients_TopCCAaverages.png');


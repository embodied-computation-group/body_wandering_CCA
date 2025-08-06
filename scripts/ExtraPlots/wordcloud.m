savepath = '/home/leah/Git/BodyWanderingCCA/figures/wordclouds/';
bw_data = readtable('/home/leah/Git/BodyWanderingCCA/data/data4correlations/BodyWanderingData.csv'); 
beh_labels = {'Future','Past','Self','Other','Pos','Neg','Words','Vivid','Vague','Spontaneous','Focus','Ruminate','Distant','Image','Arousal','Body','Breath','Heart','Movement','Bladder','Skin','Stomach'};
MDIES  = bw_data{:, beh_labels}; %
beh_labels = {'Future','Past','Self','Other','Pos','Neg','Words','Vivid','Vague','Spontaneous','Focus','Repetitive','Distant','Image','Arousal','Body','Breath','Heart','Movement','Bladder','Skin','Stomach'};
BW_av = bw_data{:, 'Body_Wandering'};  
Cog_av = bw_data{:, 'Cog_Wandering'};
Neg_diff = bw_data{:, 'Neg_Wandering'};    
BW_cor = corr(MDIES, BW_av, 'Rows', 'complete');
Cog_cor = corr(MDIES, Cog_av, 'Rows', 'complete');
Neg_cor = corr(MDIES, Neg_diff, 'Rows', 'complete');

%% Body-Wandering - MDIES items correlation word cloud %%%%%%%%%%%%%%%%%%%
% MAKE a custom colormap 
% Define your two RGB color codes
color1_rgb = [255, 112, 166]/255; % Example: light pink 
color2_rgb = [112, 214, 255]/255; % Example: light blue

% Number of colors in the colormap
num_colors = length(BW_cor); 

% Generate the custom colormap
custom_colormap = flipud([linspace(color1_rgb(1), color2_rgb(1), num_colors)', ...
                linspace(color1_rgb(2), color2_rgb(2),  num_colors)', ...
                linspace(color1_rgb(3), color2_rgb(3), num_colors)']);

% Use this custom colormap
cl = colormap(custom_colormap);

wordcloud(beh_labels', abs(BW_cor), 'Color', cl);
saveas(gcf, [savepath, 'Bodyaverage_itemCorr_WordCloud_test.png']);

%% Cog-Wandering - MDIES items correlation word cloud %%%%%%%%%%%%%%%%%%%
% MAKE a custom colormap 
% Define your two RGB color codes
color2_rgb = [112, 214, 255]/255; % Example: light blue
color1_rgb = [255, 112, 166]/255; % Example: light pink 

% Number of colors in the colormap
num_colors = length(Cog_cor); 

% Generate the custom colormap
custom_colormap = flipud([linspace(color1_rgb(1), color2_rgb(1), num_colors)', ...
                linspace(color1_rgb(2), color2_rgb(2),  num_colors)', ...
                linspace(color1_rgb(3), color2_rgb(3), num_colors)']);

% Use this custom colormap
cl = colormap(custom_colormap);

wordcloud(beh_labels', abs(Cog_cor), 'Color', cl);
saveas(gcf, [savepath, 'Cogaverage_itemCorr_WordCloud.png']);

%% Neg-Wandering - MDIES items correlation word cloud %%%%%%%%%%%%%%%%%%%
% MAKE a custom colormap 
% Define your two RGB color codes
color2_rgb = [132, 112, 255]/255; % Example: purple 
color1_rgb = [255, 214, 112]/255; % Example: light yellow

% Number of colors in the colormap
num_colors = length(Neg_cor); 

% Generate the custom colormap
custom_colormap = flipud([linspace(color1_rgb(1), color2_rgb(1), num_colors)', ...
                linspace(color1_rgb(2), color2_rgb(2),  num_colors)', ...
                linspace(color1_rgb(3), color2_rgb(3), num_colors)']);

% Use this custom colormap
cl = colormap(custom_colormap);

wordcloud(beh_labels', abs(Neg_cor), 'Color', cl);
saveas(gcf, [savepath, 'Negdiff_itemCorr_WordCloud.png']);






% edited by Leah Banellis for Gradient analysis of Body-Wandering Project 2024

% Compute designated # of permutations/spins of the input surface data
% in FreeSurfer fsaverage5.
% FORMAT SpinPermuFS(readleft,readright,permno)
% readleft     - the filename of left surface data to spin 
% readright    - the filename of right surface data to spin 
% permno       - the number of permutations
% wsname       - the name of a workspace file including all spun data to be saved
% Example   SpinPermuFS('../data/depressionFSdataL.csv','../data/depressionFSdataR.csv',100,'../data/rotationFS.mat')
% will spin prebuilt data, neurosynth map associated with 'depression', 100
% times, and save the workspace file of all spun data in ../data/rotationFS.mat
% Aaron Alexander-Bloch & Siyuan Liu 
% SpinPermuFS.m, 2018-04-22
% The implementation of generating random rotations originally described in our paper — 
% rotating the coordinates of vertices at angles uniformly chosen between zero and 360 degrees
% about each of the x (left-right), y (anterior-posterior) and z (superior-inferior) axes —
% introduces a preference towards oversampling certain rotations. 
% Thus, we modified the code to incorporate an approach, Lefèvre et al. (2018), 
% that samples uniformly from the space of possible rotations. The updated
% uniform sampling prodcedure does not require AxelRot.m anymore.
% Updated on 2018-07-18
% Update 07/31/2020 (SMW): will automatically remove medial wall for
% fsaverage5. may need to change if not fsaverage5 (10242 vertices per
% hemisphere)
% edited by Leah Banellis for Gradient analysis of Body-Wandering Project 2024

%% initialise variables and paths
clear all
%Set up paths
fshome = '/Applications/freesurfer/7.4.1/'; %getenv('FREESURFER_HOME');
fsmatlab = sprintf('%s/matlab',fshome);
path(path,fsmatlab);
addpath('/Users/au704655/Documents/Body_wandering/Scripts/spin-test-master/scripts')

permno = 2500; %permutation starts
%gradientno = 5; %number gradient atlases

%% load gradient atlas and parcellation atlas in surface space
% schaefer-subcort combined atlas in surface space
lh_parcels = load_mgh('/Users/au704655/Documents/Body_wandering/Data/gradients/schaefer200subcort_atlas.L.mgh');
rh_parcels = load_mgh('/Users/au704655/Documents/Body_wandering/Data/gradients/schaefer200subcort_atlas.R.mgh');
% Combine left and right hemisphere parcels
surface_parcels = [lh_parcels; rh_parcels];

% gradient atlas in surface space
datal = load_mgh('/Users/au704655/Documents/Body_wandering/Data/gradients/gradientcort_surface/corticalsubcortical/gradient1/gradient1cort_surface.L.mgh');
datar =  load_mgh('/Users/au704655/Documents/Body_wandering/Data/gradients/gradientcort_surface/corticalsubcortical/gradient1/gradient1cort_surface.R.mgh');
% Combine left and right gradient atlas
gradients_data = [datal; datar];

%% Calculate average across gradient values for each parcel (but assess top 1% of CCA parcels)
unique_parcels = unique(surface_parcels); 
for i = 1:length(unique_parcels)
    parcel = unique_parcels(i);
    vertex_indices = find(surface_parcels == parcel);
    average_gradients(i,1) = parcel; %parcel num
    average_gradients(i,2) = mean(gradients_data(vertex_indices)); % average gradient for that parcel
end
% % check surface-space averages are similar to parcel averages in python
% csv_data = readmatrix('/Users/au704655/Documents/Body_wandering/Data/gradients/gradient_averages_cca.csv'); %readtable
% % Match the entries and extract corresponding second column values
% [common_vals, csv_indices, gradient_indices] = intersect(csv_data(:,1), average_gradients(:,1));
% % Perform the correlation to check 
% corr(csv_data(csv_indices,2), average_gradients(gradient_indices,2))
% scatter(csv_data(csv_indices,2), average_gradients(gradient_indices,2))

%% load freesurface files for spin-test
%extract the vertice coordinates from the corresponding surface for plotting before/after spin-test
[verticespl, facespl] = freesurfer_read_surf(fullfile(fshome,'subjects/fsaverage/surf/lh.pial'));
[verticespr, facespr] = freesurfer_read_surf(fullfile(fshome,'subjects/fsaverage/surf/rh.pial'));
% plot graident atlas before spin-test
custommap=colormap('jet');
subplot(2,2,1);
mincol=min(datal);
maxcol=max(datal);
plotFSsurf(facespl,verticespl,datal,custommap,mincol,maxcol,[-90 0]);
title('Lateral View of Initial Left');
subplot(2,2,2);
plotFSsurf(facespl,verticespl,datal,custommap,mincol,maxcol,[90 0]);
title('Medial View of Initial Left');

% left % right annotation files:
%[vl, left_labels, ctl] = read_annotation(fullfile(fshome,'/subjects/fsaverage/label/lh.aparc.a2009s.annot'));
%[vr,right_labels,ctr] = read_annotation(fullfile(fshome,'/subjects/fsaverage/label/rh.aparc.a2009s.annot'));

%%extract the corresponding sphere surface coordinates for rotation
[verticesl, ~] = freesurfer_read_surf(fullfile(fshome,'subjects/fsaverage/surf/lh.sphere'));
[verticesr, ~] = freesurfer_read_surf(fullfile(fshome,'subjects/fsaverage/surf/rh.sphere'));

%% Spin-Test
%Use rng to initialize the random generator for reproducible results.
rng(0);
%initialize variables to save rotation
bigrotl=[];
bigrotr=[];
%distfun = @(a,b) sqrt(bsxfun(@minus,bsxfun(@plus,sum(a.^2,2),sum(b.^2,1)),2*(a*b)));
%function to calculate Euclidian distance, deprecated 2019-06-18 see home page
I1 = eye(3,3);
I1(1,1)=-1;
bl=verticesl;
br=verticesr;

% run spin-test for number of permutations
average_gradients_perm = zeros(length(unique(surface_parcels)), permno); % ..gradientno,
for j=1:permno
    j
    %the updated uniform sampling procedure
    A = normrnd(0,1,3,3);
    [TL, temp] = qr(A);
    TL = TL * diag(sign(diag(temp)));
    if(det(TL)<0)
        TL(:,1) = -TL(:,1);
    end
    %reflect across the Y-Z plane for right hemisphere
    TR = I1 * TL * I1;
    bl =bl*TL;
    br = br*TR;    
    
    %Find the pair of matched vertices with the min distance and reassign
    %values to the rotated surface.
    %distl=distfun(verticesl,bl'); % deprecated 2019-06-18 see home page
    %distr=distfun(verticesr,br'); % deprecated 2019-06-18 see home page
    %[~, Il]=min(distl,[],2); % deprecated 2019-06-18 see home page
    %[~, Ir]=min(distr,[],2); % deprecated 2019-06-18 see home page
    Il = nearestneighbour(verticesl', bl'); % added 2019-06-18 see home page
    Ir = nearestneighbour(verticesr', br'); % added 2019-06-18 see home page

    %save rotated data
    bigrotl= datal(Il)'; %[bigrotl; datal(Il)'];
    bigrotr= datar(Ir)'; %[bigrotr; datar(Ir)'];
    % it is also feasible to save Il Ir and apply them to different datasets
    % for repeated use
    %If annotation file is used, annotation file for each rotation could be
    %saved by write_annotation.m of FreeSurfer

%     % plot gradient altas after spin-test
%     %figure;
%     subplot(2,2,3);
%     mincol=min(bigrotl);
%     maxcol=max(bigrotl);
%     plotFSsurf(facespl,verticespl,bigrotl,custommap,mincol,maxcol,[-90 0]);
%     title('Lateral View of resampled Left');
%     subplot(2,2,4);
%     plotFSsurf(facespl,verticespl,bigrotl,custommap,mincol,maxcol,[90 0]);
%     title('Medial View of resampled Left');

    %% Store 'spinned' gradient averages for each spin (for null-distribution)
    % Combine rotated left and right gradient data
    gradients_data_rot = [bigrotl'; bigrotr'];
    % Calculate average gradient values for each parcel
    for i = 1:length(unique_parcels)
        parcel = unique_parcels(i);
        vertex_indices = find(surface_parcels == parcel);
        average_gradients_perm(i,j) = mean(gradients_data_rot(vertex_indices));
    end
end
gradient_spin.average_gradients_perm = average_gradients_perm;
gradient_spin.parcelnums = unique_parcels;
gradient_spin.average_gradients = average_gradients;

% save gradient averages and spin-test gradient averages
save('gradient1_spin.mat', 'gradient_spin');

%save(wsname,'bigrotl','bigrotr')
%save bigrotl and bigrotr in a workspace file for the null distribution
%use it in pvalvsNull.m to caclulate pvalue

%% plots QC %%
% % check surface coversion of schaefer atlas worked:
% lh_parcels = load_mgh('/Users/au704655/Documents/Body_wandering/Data/gradients/schaefer200subcort_atlas.L.mgh');
% rh_parcels = load_mgh('/Users/au704655/Documents/Body_wandering/Data/gradients/schaefer200subcort_atlas.R.mgh');
% % change parcel numbers so scale for plotting is better
% mapping = containers.Map(unique(lh_parcels), 1:length(unique(lh_parcels)));
% lh_parcels_sequential = arrayfun(@(x) mapping(x), lh_parcels);
% mapping = containers.Map(unique(rh_parcels), 1:length(unique(rh_parcels)));
% rh_parcels_sequential = arrayfun(@(x) mapping(x), rh_parcels);
% 
% custommap=colormap('jet');
% subplot(2,2,1);
% mincol=1; %min(lh_labels); (subcort vars (3000-7007) much higher parcel nums than cortical)
% maxcol=107; %max(lh_labels);
% plotFSsurf(facespl,verticespl,lh_parcels_sequential,custommap,mincol,maxcol,[-90 0]);
% title('Schaefer atlas Left');
% subplot(2,2,2);
% plotFSsurf(facespl,verticespl,lh_parcels_sequential,custommap,mincol,maxcol,[90 0]);
% title('Schaefer atlas Left');
% 
% subplot(2,2,3);
% mincol=1;
% maxcol=108;
% plotFSsurf(facespl,verticespl,rh_parcels_sequential,custommap,mincol,maxcol,[-90 0]);
% title('Schaefer atlas Right');
% subplot(2,2,4);
% plotFSsurf(facespl,verticespl,rh_parcels_sequential,custommap,mincol,maxcol,[90 0]);
% title('Schaefer atlas Right');

 
% % create mask of ROI (1003 = Thalamus L in this example)
% lh_mask = zeros(size(lh_labels));
% lh_mask(lh_labels == 6000) = 1000; % 1000 just so a very different colour to 0
% 
% custommap=colormap('jet');
% subplot(2,2,1);
% mincol=min(datal);
% maxcol=max(datal);
% plotFSsurf(facespl,verticespl,lh_mask,custommap,mincol,maxcol,[-90 0]);
% title('Thalamus Left');
% subplot(2,2,2);
% plotFSsurf(facespl,verticespl,lh_mask,custommap,mincol,maxcol,[90 0]);
% title('Thalamus Left');
% 
% % create mask of ROI (1003 = Lateral Occipital Cortex Inf L in this example)
% lh_mask = zeros(size(lh_labels));
% lh_mask(lh_labels == 1003) = 1000; % 1000 just so a very different colour to 0
% 
% subplot(2,2,3);
% mincol=min(datal);
% maxcol=max(datal);
% plotFSsurf(facespl,verticespl,lh_mask,custommap,mincol,maxcol,[-90 0]);
% title('Lateral Occipital Cortex Inf L');
% subplot(2,2,4);
% plotFSsurf(facespl,verticespl,lh_mask,custommap,mincol,maxcol,[90 0]);
% title('Lateral Occipital Cortex Inf L');


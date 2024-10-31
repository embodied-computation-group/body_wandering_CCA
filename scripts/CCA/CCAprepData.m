function [X1,Y1,Z1] = A_CCAprepData_BrainAbend_BWitems()
% prepare netmats (vectorised connectivity data), subject data (body-wandering), and nuisance vars for CCA

%   Outputs
% X1  sub variables of interest (i.e., Mind/Body Wandering Items)
% Y1  netmats (Connectivity Data)
% Z1  confounds 

addpath(genpath("/Users/au704655/Documents/Body_wandering/Scripts/CCA/supporting_packages/"))

%% Key flags     
saveOutput          = 1;

%% Data directories
dirImaging         	= '/Users/au704655/Documents/Body_wandering/Data/abend_schaefer200_16subcort/correlation'; % connectivity data
dirTask             = '/Users/au704655/Documents/Body_wandering/Data'; % task data location
subLevelFile        = 'mw_rest_probes.csv';                 % task data filename (mind wandering items)
dirConfounds        = '/Users/au704655/Documents/Body_wandering/Data/demographics.csv';  % nuisance data location

%% ------------------------------------------------
% Load FUNCTIONAL CONNECTIVITY data and create NetMats matrix
% -------------------------------------------------
subjectFCfiles = dir(fullfile(dirImaging, '*.npy'));
netMats             = [];

for thisFile = 1 : length(subjectFCfiles)  
    % Load .npy matrices
    subFilePath     = fullfile(dirImaging, subjectFCfiles(thisFile).name);
    thisSubMat      = readNPY(subFilePath);
    
    % Reshape to row vector
    thisSubArray    = upperMatTri2Vector(thisSubMat);
    
    % Get subject ID and add to first column
    thisSubID       = str2double(subFilePath(end-7:end-4));
    thisSubArray    = [thisSubID, thisSubArray];
    
    % Add to multi-sub matrix
    netMats = [netMats; thisSubArray];
end
% netMats (subjects, long array of FC values)

%% -----------------------------------------------
% Load TASK subject-level data - mind/body wandering factors
% ------------------------------------------------
dataPathTask        = fullfile(dirTask,subLevelFile);
taskDataTable       = readtable(dataPathTask,'TreatAsEmpty',{'NA'}); %struct2table(tdfread(dataPathTask)); 
taskDataTable(:,"Var1")  = []; % remove index col

% change participant_id from string to numeric
% load subject IDs for body-wandering data
subIDs = readtable('/Users/au704655/Documents/Body_wandering/Data/subjectIDs_4_EFA_oblimin_PA_loadings.csv');
subIDs= char(table2array(subIDs(:,2)));
subIDs_num=[];
for n = 1:size(subIDs,1)
    subIDs_num(n) = str2double(subIDs(n,5:8));
end

taskDataTable = [array2table(subIDs_num'),taskDataTable];
taskDataTable = renamevars(taskDataTable,"Var1","participant_id");

% remove spontaneous - low KMO
taskDataTable(:,"Spontaneous") = [];

taskDataArray       = table2array(taskDataTable);

% remove subjects with NaN's
taskDataArray(any(isnan(taskDataArray), 2), :) = [];

%% Clean up input matrices

% Remove subjects who don't have both data types
usableSubs          = intersect(netMats(:,1),taskDataArray(:,1)); %intersect(netMats(:,1),taskDataArray(:,2));

% Get only usable subjects from netMats and task 
Y0                  = netMats(ismember(netMats(:,1),usableSubs),:);
X0                  = taskDataArray(ismember(taskDataArray(:,1),usableSubs),:); %taskDataArray(ismember(taskDataArray(:,2),usableSubs),:);

Y0 = sortrows(Y0,1);
X0 = sortrows(X0,1);

%% ------------------------------------------------
% Load CONFOUNDS
% -------------------------------------------------
confoundsDataTable    	= readtable(dirConfounds);      

% extract confound variables of interest
ID              = cellfun(@(x){x(5:8)}, confoundsDataTable.participant_id);
gender          = confoundsDataTable.gender;              % cannot have binary variable
age             = confoundsDataTable.age;
weight          = confoundsDataTable.weight;
height          = confoundsDataTable.height;

% calculate bmi (height in cm, weight in kilos)(BMI = kg/m2)
bmi = (weight./(height.*height))*10000;

% define vmp1/2 session for confound array
session = ones(length(ID),1);
session(str2num(cell2mat(ID)) < 263) = -1;

confoundsArray 	= horzcat(str2num(cell2mat(ID)),age,bmi,gender,session);%
% remove gender = 3 (other)
confoundsArray(find(confoundsArray(:,4) == 3),4) = NaN;
% make gender 1 & -1
confoundsArray(find(confoundsArray(:,4) == 2),4) = -1;

% match subjects
confoundsArray  = sortrows(confoundsArray,1);
Z0              = confoundsArray(ismember(confoundsArray(:,1),usableSubs),:);  

% Find subjects without phenotype data and remove from ccaReadySubVars &
% ccaReadyNetMats
[row, col] = find(isnan(Z0)); %setdiff(X0(:,1),Z0(:,1));
missingPhenoSubs    = Z0(row,1);
X0                  = X0(~ismember(X0(:,1),missingPhenoSubs),:);
Y0                  = Y0(~ismember(Y0(:,1),missingPhenoSubs),:);
Z0(row,:) = [];

%% Remove subject IDs
X1      = X0(:,2:end);
Y1      = Y0(:,2:end);
Z1      = Z0(:,2:end);

% convert Y1 array to double
Y1 = double(Y1);

%% Save the input matrices
if saveOutput
    filename    = 'cca_inputs_abend216_BodywanderingItems';
    dirCCAinput    	= fullfile('.','input',filename);

    %% for CCA/PLS toolbox
    X = Y1; % connectivity
    Y = X1; % beh data
    C = Z1; % confounds
    save(dirCCAinput, 'X','Y','C')
end


%% Plot CCA inputs
savepath = '/Users/au704655/Documents/Body_wandering/Figures/CCA/Abend216_Items/CCAinput/'; 
labels = {'Future', 'Past', 'Self', 'Other', 'Pos', 'Neg', 'Words', 'Vivid', 'Vague', 'Focus', 'Ruminate', 'Distant', 'Image', 'Arousal', 'Body', 'Breath', 'Heart', 'Movement', 'Bladder', 'Skin', 'Stomach'};  

figure; Xplot = imagesc(X);
saveas(Xplot, [savepath, 'X_heatmap.png'], 'png');
figure; Yplot = imagesc(Y);
saveas(Yplot, [savepath, 'Y_heatmap.png'], 'png');
figure; Cplot = imagesc(C);
saveas(Cplot, [savepath, 'C_heatmap.png'], 'png');

%% check correlations
[R, p_values] = corrcoef(Y);

num_labels = numel(labels);
figure('Position', [100, 100, ceil(num_labels/4) * 150, ceil(num_labels/4) * 150]);
% colour by significance < 0.05
corrmatrix = imagesc(R, 'AlphaData', p_values < 0.05);
colorbar;
title('Correlation Matrix');
% r-values on top
[row, col] = size(R);
for i = 1:row
    for j = 1:col
        text(j, i, num2str(R(i, j), '%.2f'), 'HorizontalAlignment', 'center');
    end
end
% labels
xticks(1:length(R));
yticks(1:length(R));
xticklabels(labels);
yticklabels(labels);
% colours
colormap('cool')
caxis([-1, 1]);
% Save the image as PNG
saveas(corrmatrix, [savepath, 'Beh_corrmatrix.png'], 'png');

%% save histograms of Y in 1 fig
figure('Position', [100, 100, ceil(num_labels/4) * 150, ceil(num_labels/4) * 150]);
num_cols = size(Y, 2);
num_rows = ceil(num_cols / 4);
for n = 1:num_cols
    subplot(num_rows, ceil(num_cols / num_rows), n);
    hist(Y(:, n));
    xlabel(labels(n));
    ylabel('Frequency');
end
sgtitle('Histograms of Variables in Y');
saveas(gcf, [savepath, 'Y_histograms.png'], 'png');


end



close all; clear;
mice = {
    %'5038_2L'; Plastic bar
    %'5038_2L1R'; Plasric bar
    %'5038_2L2R'; Plastic bar
    %'5039_1L'; Plastic bar
    %'5039_1L1R'; Plastic bar
    %'1442_NM'; Plastic bar
    
    '0931_1B','NH',1;%No holder
    '0931_1L','NH',1;%No holder
    '0931_1R','NH',2; %No holder
    '1675_1R','NH',2;%No holder
    
    '1442_1L1R','H',2;%Holder
    '9145_1L','H',2;%Holder
    'K81_1R','H',1; % Holder
    'K81_1L','H',2; % Holder
    'K81_1L1R','H',2; % Holder
    };

miceTumors = {'NH','NH','NH','NH','NH','NH','H','H','H','H','H','H','H','H','H'}
miceTumorsId = [1,1,1,1,1,1,0,0,0,0,0,0,0,0,0];


diceCoefficientsForBodyOverlay = [];
diceCoefficientsForTumorRegistration = [];
bodyAreas = [];

correlationDCEandOEOT = []

for k=1:length(mice)
    load(['/media/gehrun01/Dropbox1/Dropbox/Cloud/CRUK CI/Masters Thesis/Framework/thesis-db/mat/' mice{k,1} '.mat']);
    t2anatomical = [];
    dceMri = [];
    if(size(mouse.mri.t2,3)>1)
        t2anatomical = flip(mouse.mri.t2(:,:,2),2);
    else
        t2anatomical = flip(mouse.mri.t2,2);
    end
    
    if isfield(mouse.mri,'dce') && k ~= 9
        size(mouse.mri.dce)
        if(size(mouse.mri.dce,3)>1)
            dceMri = flip(mouse.mri.dce(:,:,2),2);
        else
            dceMri = flip(mouse.mri.dce,2);
        end 
        dceMri = imresize(dceMri,size(t2anatomical));
    end
    
    
    
    
     mouseMsot = mouse.msp.totalHbUnderAir;
     mouseMsotOE = mouse.msp.oxygenSaturationDelta;
     if k == 5 || k == 6
         mouseMsot = flip(mouseMsot,2);
         mouseMsotOE = flip(mouseMsotOE,2);
     end
     

    
    load(['resources/body-rois/' mice{k,1} '.mat']);
    % Resolution correction
    % MSOT
    % 267x267 @ 20mm FOV: 75 micron/pixel
    % 334x334 @ 25mm FOV: 75 micron/pixel
    % MRI
    % 256x256 @ 40mm FOV: 156 micron/pixel
    %
    % Scaling factor for MSOT is 75/156=0.4808
    scalingFactor = 0.4808;
    if k==5
       %figure;imshow(mouseMsot,[])
      % figure;imshow(t2anatomical,[])
    end
    mouseMsot = imresize(mouseMsot,0.4808);
    %mouseMsotOE = imresize(mouseMsotOE,0.4808);
    bodyRoi.msot = imresize(bodyRoi.msot,size(mouseMsot));
    bodyRoi.mri = imresize(bodyRoi.mri,size(t2anatomical));
    %t2anatomical = imresize(t2anatomical,0.78);
  
    %figure;imshowpair(mouseMsot,t2anatomical,'montage')

    load(['resources/registration-landmarks/' mice{k,1} '.mat']);
    
    
    if size(mouseMsot)==[129 129]
    landmarkPoints.msot = landmarkPoints.msot.*1.29;
    end
    
    if size(mouseMsot)==[161 161]
    landmarkPoints.msot = landmarkPoints.msot.*1.28;
    end
    
    landmarkPoints.mri = landmarkPoints.mri.*1.28;
    
    
    
    %%figure;showMatchedFeatures(mouseMsot,t2anatomical,landmarkPoints.msot,landmarkPoints.mri,'montage');
    
    pointsetMriMsotRegistration = fitgeotrans(landmarkPoints.mri, landmarkPoints.msot, 'similarity');
    
    
    
    pointsetMriRegistered = imwarp(t2anatomical,pointsetMriMsotRegistration,'OutputView',imref2d(size(mouseMsot)));
     if isfield(mouse.mri,'dce') && k ~= 9
    pointsetMriDCERegistered = imwarp(dceMri,pointsetMriMsotRegistration,'OutputView',imref2d(size(mouseMsot)));
     end
    pointsetMriBinRegistered = imwarp(bodyRoi.mri,pointsetMriMsotRegistration,'OutputView',imref2d(size(bodyRoi.msot)));
    %figure;imshowpair(pointsetMriBinRegistered,bodyRoi.msot);
    diceCoefficientsForBodyOverlay(k,:) = [calculateDiceSimilarityCoefficient(bodyRoi.mri,bodyRoi.msot),calculateDiceSimilarityCoefficient(pointsetMriBinRegistered,bodyRoi.msot)];
    bodyAreas(k,:) = [nnz(pointsetMriBinRegistered),nnz(bodyRoi.msot)];
    
    %figure;imshowpair(mat2gray(pointsetMriRegistered),mat2gray(mouseMsot));
    
    for l=1:mice{k,3}
        
        load(['resources/tumor-registration-landmarks/tumor_landmarks_' mice{k,1} '_' num2str(l) '.mat'])
        load(['resources/tumor-registration-contours/' mice{k,1} '_msot_' num2str(l) '.mat']);
        tumorRoiMsot = tumorRoi;
        load(['resources/tumor-registration-contours/' mice{k,1} '_mri_' num2str(l) '.mat']);
        tumorRoiMri = tumorRoi;
        %figure;imshow(tumorRoiMri)
        %figure;imshow(tumorRoiMsot)
        pointsetPostMriMsotRegistration = fitgeotrans(tumorLandmarks.mri, tumorLandmarks.msot, 'similarity');
        pointsetTumorMriBinRegistered = imwarp(tumorRoiMri,pointsetPostMriMsotRegistration,'OutputView',imref2d(size(tumorRoiMsot)));
        pointsetTumorMriRegistered = imwarp(pointsetMriRegistered,pointsetPostMriMsotRegistration,'OutputView',imref2d(size(mouseMsot)));
         if isfield(mouse.mri,'dce') && k ~= 9
        pointsetTumorMriDCERegistered = imwarp(pointsetMriDCERegistered,pointsetPostMriMsotRegistration,'OutputView',imref2d(size(mouseMsot)));
         end
        diceCoefficientsForTumorRegistration(end+1,:) = [calculateDiceSimilarityCoefficient(tumorRoiMri,tumorRoiMsot),calculateDiceSimilarityCoefficient(pointsetTumorMriBinRegistered,tumorRoiMsot)];
        
        if k>6 && k<10
            %figure;imshow(pointsetTumorMriRegistered,[])
            %figure;imshow(mouseMsot,[])
            %figure;imshowpair(mat2gray(pointsetTumorMriRegistered),mat2gray(mouseMsot));title(k);
        end
        
        %figure;imshowpair(mat2gray(pointsetTumorMriRegistered),mat2gray(mouseMsot),'montage');title(k);
       % figure;showMatchedFeatures(tumorRoiMri,tumorRoiMsot,tumorLandmarks.mri, tumorLandmarks.msot,'montage')
        % figure;imshowpair(tumorRoiMri,tumorRoiMsot);
        %figure;imshowpair(pointsetTumorMriBinRegistered,tumorRoiMsot);
        
        if isfield(mouse.mri,'dce') && k ~= 9
            disp('DCE exists')
            figure;subplot(1,2,1);imshow(pointsetTumorMriDCERegistered);
            subplot(1,2,2);imshow(mouseMsotOE);
            correlationDCEandOEOT(end+1) = corr(pointsetTumorMriDCERegistered(tumorRoiMsot),mouseMsotOE(tumorRoiMsot),'Rows','complete')
            %correlationDCEandOEOT(end)
        end
        
        
    end
    
    
    
end

%% Before/after rough registration
figure;boxplot(diceCoefficientsForBodyOverlay(:,1),mice(:,2),'Colors',[0 0 0;0 0 0])
%title('No body alignment')
set(findobj(gca,'type','line'),'linew',2)
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');
plt = Plot();
plt.ShowBox = false;
plt.BoxDim = [7 7];
plt.YLim = [0.4 1];
plt.XMinorTick = false;
plt.TickLength = [0.01 0.01];
plt.export('noBodyReg.pdf');

%Print stats
%unpaired t test
median(diceCoefficientsForBodyOverlay(1:4,1))
iqr(diceCoefficientsForBodyOverlay(1:4,1))
median(diceCoefficientsForBodyOverlay(5:9,1))
iqr(diceCoefficientsForBodyOverlay(5:9,1))
[h,p]=ttest2(diceCoefficientsForBodyOverlay(1:4,1),diceCoefficientsForBodyOverlay(5:9,1))


%%
figure;boxplot(diceCoefficientsForBodyOverlay(:,2),mice(:,2),'Colors',[0 0 0;0 0 0])
set(findobj(gca,'type','line'),'linew',2)
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');
plt = Plot();
plt.ShowBox = false;
plt.BoxDim = [7 7];
plt.YLim = [0.4 1];
plt.XMinorTick = false;
plt.TickLength = [0.01 0.01];
plt.export('withBodyReg.pdf');

%Print stats
%unpaired t test
median(diceCoefficientsForBodyOverlay(1:4,2))
iqr(diceCoefficientsForBodyOverlay(1:4,2))
median(diceCoefficientsForBodyOverlay(5:9,2))
iqr(diceCoefficientsForBodyOverlay(5:9,2))
[h,p]=ttest2(diceCoefficientsForBodyOverlay(1:4,2),diceCoefficientsForBodyOverlay(5:9,2))


%% Between tests
[h,p]=ttest(diceCoefficientsForBodyOverlay(1:4,1),diceCoefficientsForBodyOverlay(1:4,2))

[h,p]=ttest(diceCoefficientsForBodyOverlay(5:9,1),diceCoefficientsForBodyOverlay(5:9,2))



%% Before after software refinement
figure;boxplot(diceCoefficientsForTumorRegistration,{'Before','After'},'Colors',[0 0 0;0 0 0])
set(findobj(gca,'type','line'),'linew',2)
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');
plt = Plot();
plt.ShowBox = false;
plt.BoxDim = [7 7];
plt.YLim = [0.4 1];
plt.XMinorTick = false;
plt.TickLength = [0.01 0.01]
plt.export('softwareReg.pdf');

% paired t-test
median(diceCoefficientsForTumorRegistration(:,1))
iqr(diceCoefficientsForTumorRegistration(:,1))
median(diceCoefficientsForTumorRegistration(:,2))
iqr(diceCoefficientsForTumorRegistration(:,2))
[h,p]=ttest(diceCoefficientsForTumorRegistration(:,1),diceCoefficientsForTumorRegistration(:,2))

%% T2 star 

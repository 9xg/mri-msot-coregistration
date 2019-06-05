close all

otMriPairs = {};
mouseMriOt = [];
mriOTTumors = {
    '5038_2L','',2,'GEDA','Control';
    '5038_2L1R','',1,'GEDA','Control';
    '5038_2L2R','',1,'GEDA','Control';
    '5039_1L','',1,'GEDA','Control';
    '5039_1L1R','',2,'GEDA','Control';
    '1442_NM','',1,'PC3','Atovaquone_Control';
    '1442_1L1R','',1,'PC3','Atovaquone_Control';
    '9145_1L','',1,'PC3','Atovaquone_Treated';
    };

for k=1:length(mriOTTumors)
    load(['/media/gehrun01/Dropbox1/Dropbox/Cloud/CRUK CI/Masters Thesis/Framework/thesis-db/mat/' mriOTTumors{k,1} '.mat']);
    otMriPairs{end+1,1} = mouse.msp.totalHbUnderAir;
    if(size(mouse.mri.t2,3)>1)
        otMriPairs{end,2} = flip(imresize(mouse.mri.t2(:,:,2),0.5),2);
    else
        otMriPairs{end,2} = flip(mouse.mri.t2,2);
    end
    otMriPairs{end,3} = mriOTTumors{k};
    
    mouse.msp = rmfield(mouse.msp,'Air');
    mouse.msp = rmfield(mouse.msp,'O2');
    
    mouseMriOt(end+1).ot = mouse.msp;
    mouseMriOt(end).ot.icg = mouse.icg.icgAmp;
    mouseMriOt(end).maska = mouse.tumor;
    mouseMriOt(end).mri = mouse.mri;

end
%%
mriOTTumors2 = {
    '0931_1B';
    '0931_1L';
    '0931_1R';
    '1675_1R';
    };
for k=1:length(mriOTTumors2)
    load(['/media/gehrun01/Dropbox1/Dropbox/Cloud/CRUK CI/Masters Thesis/Framework/thesis-db/old-mri-msot-registration/MRI/' mriOTTumors2{k,1} '.mat']);
    load(['/media/gehrun01/Dropbox1/Dropbox/Cloud/CRUK CI/Masters Thesis/Framework/thesis-db/old-mri-msot-registration/OT/' mriOTTumors2{k,1} '.mat']);
    otMriPairs{end+1,1} = MSP_Air(:,:,1)+MSP_Air(:,:,2);
    otMriPairs{end,2} = flip(squeeze(T2_data(:,:,2)),2);
    otMriPairs{end,3} = mriOTTumors2{k};

end


%%
pointsetMriMsotMetrics = [];
for o=1:length(otMriPairs)
    %figure;imshow(otMriPairs{o,2},[]);
    load(['resources/msot-mri-body-contours/' otMriPairs{o,3} '.mat'],'bodyRoiMask');
    mriMask = imfill(imbinarize(mat2gray(otMriPairs{o,2}),0.1),'holes');
    msotMask = bodyRoiMask;
    %figure;imshow(mriMask)
    %figure;imshow(msotMask)
    msotMask = bodyRoiMask;
    if(o>9 || o<6)
        msotMask = imresize(msotMask,0.5);
        otMriPairs{o,1} = imresize(otMriPairs{o,1},0.5);
    end
    
    
    mriMetric = otMriPairs{o,2};
    % Point set registration
    load(['/media/gehrun01/Dropbox1/Dropbox/Cloud/CRUK CI/Masters Thesis/Framework/thesis-db/coregistration-parameters/mri-ot-per-mouse-registration/' otMriPairs{o,3} '_cP.mat']);
    if(size(coregistrationParameters.mriPoints,1)>2)
        pointsetMriMsotRegistration = fitgeotrans(coregistrationParameters.mriPoints, coregistrationParameters.msotPoints, 'similarity');
    else
        pointsetMriMsotRegistration = fitgeotrans(coregistrationParameters.mriPoints, coregistrationParameters.msotPoints, 'nonreflectivesimilarity');
    end
% 
%     if o==6
%     figure;imshowpair(mriMask,otMriPairs{o,1})
%     figure;imshow(mriMask)
%     figure;imshow(msotMask)
%     figure;showMatchedFeatures(mriMask,msotMask,coregistrationParameters.mriPoints,coregistrationParameters.msotPoints);
%    end
    pointsetMriBinRegistered = imwarp(mriMask,pointsetMriMsotRegistration,'OutputView',imref2d(size(otMriPairs{o,1})));
    pointsetMriRegistered = imwarp(mriMetric,pointsetMriMsotRegistration,'OutputView',imref2d(size(otMriPairs{o,1})));
    pointsetMriMsotMetrics(o,:) = [calculateDiceSimilarityCoefficient(pointsetMriBinRegistered,msotMask),immse(double(pointsetMriRegistered),double(otMriPairs{o,1})),nnz(pointsetMriBinRegistered)/nnz(msotMask)];

  figure;imshowpair(pointsetMriBinRegistered,msotMask);

end

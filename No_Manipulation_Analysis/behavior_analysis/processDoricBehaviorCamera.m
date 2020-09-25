function processDoricBehaviorCamera(VideoDirectory,roiDir,folderTracking)

%% Create ROI(s) for each imaging session
defineBehaviorROIs(VideoDirectory,roiDir);
%% Track Movement withing each ROI
Track_BehaviorCam(VideoFolder,roiDIr,folderTracking);
%% Binarize Behaviors
cd(folderTracking);
trackFiles=dir('*behaviorTracking.mat');
for filNum=1:size(trackFiles,1)
    saveFlag='n';
    load(trackFiles(filNum).name);
    roiNames=fieldnames(trackCam.Tracking);
    for roiNum=1:size(roiNames,1)
        restDur=5*trackCam.Params.frameRate(1);
        eventDur=2*trackCam.Params.frameRate(1);
        if ~isfield(trackCam,'binarizedData')
            [binBehavior,new_eventInds,movementThresh]=binarize_BehaviorCam(trackCam.Tracking.(roiName{roiNum}).smoothTracking,restDur,eventDur);
            trackCam.binarizedData.(roiName{roiNum}).binBehavior=binBehavior;
            trackCam.binarizedData.(roiName{roiNum}).behaviorInds=new_eventInds;
            trackCam.binarizedData.(roiName{roiNum}).behaviorThresh=movementThresh;
            saveFlag='y';
        else
            theFlag=input('Do you want to reanalyze camera data binarization? (y/n)','s');
            if strcmpi(theFlag,'y')
                [binBehavior,new_eventInds,movementThresh]=binarize_BehaviorCam(trackCam.Tracking.(roiName{roiNum}).smoothTracking,restDur,eventDur);
                trackCam.binarizedData.(roiName{roiNum}).binBehavior=binBehavior;
                trackCam.binarizedData.(roiName{roiNum}).behaviorInds=new_eventInds;
                trackCam.binarizedData.(roiName{roiNum}).behaviorThresh=movementThresh;
                saveFlag='y';
            end
        end
    end
    if strcmpi(saveFlag,'y')
    save(trackFiles(filNum).name,'trackCam','-v7.3');
    end
end

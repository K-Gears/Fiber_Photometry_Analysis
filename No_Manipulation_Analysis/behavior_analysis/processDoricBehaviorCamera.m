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
    if ~isfield(trackCam,'binarizedData')
        for roiNum=1:size(roiNames,1)
            restDur=5*trackCam.Params.frameRate(1);
            eventDur=2*trackCam.Params.frameRate(1);
            roiName=roiNames{roiNum};
            [binBehavior,new_eventInds,movementThresh]=binarize_BehaviorCam(trackCam.Tracking.(roiNames{roiNum}).smoothTracking,restDur,eventDur,roiName);
            trackCam.binarizedData.(roiNames{roiNum}).binBehavior=binBehavior;
            trackCam.binarizedData.(roiNames{roiNum}).behaviorInds=new_eventInds;
            trackCam.binarizedData.(roiNames{roiNum}).behaviorThresh=movementThresh;
            saveFlag='y';
        end
    else
        theFlag=input('Do you want to reanalyze camera data binarization? (y/n)','s');
        if strcmpi(theFlag,'y')
            for roiNum=1:size(roiNames,1)
                restDur=5*trackCam.Params.frameRate(1);
                eventDur=2*trackCam.Params.frameRate(1);
                roiName=roiNames{roiNum};
                [binBehavior,new_eventInds,movementThresh]=binarize_BehaviorCam(trackCam.Tracking.(roiNames{roiNum}).smoothTracking,restDur,eventDur,roiName);
                trackCam.binarizedData.(roiNames{roiNum}).binBehavior=binBehavior;
                trackCam.binarizedData.(roiNames{roiNum}).behaviorInds=new_eventInds;
                trackCam.binarizedData.(roiNames{roiNum}).behaviorThresh=movementThresh;
                saveFlag='y';
            end
        end
    end
    if strcmpi(saveFlag,'y')
        save(trackFiles(filNum).name,'trackCam','-v7.3');
    end
end
   
end

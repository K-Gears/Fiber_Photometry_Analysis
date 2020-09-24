function Track_BehaviorCam(VideoFolder,folderROIs,folderTracking)
% This function loads in .avi files collected using BFLY usb cameras,
% and a user defined ROI generated in 'defineBehaviorROIs.m' and calculates
% the sum of the absolute difference in pixel intensity within the ROI as a proxy for
% movement.
% INPUTS:
%VideoFolder: directory containing behavior camera avi files
%folderROIs: directory containing user defined ROIs for each imaging
%session
%folderTracking: directory to save structure 'trackCam' containing movement
%for each ROI to be binarized in ExtractFPData.

%% Find behavior camera videos
cd(VideoFolder);
subfolders=dir;
subfolders(~[subfolders.isdir])=[];
tf=ismember({subfolders.name},{'.','..'});
subfolders(tf)=[];
for folderNum=1:size(subfolders,1)
    cd([subfolders(folderNum).folder '\' subfolders(folderNum).name]);
    anFolders=dir;
    anFolders(~[anFolders.isdir])=[];
    tf=ismember({anFolders.name},{'.','..'});
    anFolders(tf)=[];
    for anNum=1:size(anFolders,1)
        close all
        cd([anFolders(anNum).folder '\' anFolders(anNum).name]);
        thebreaks=strfind(anFolders(anNum).folder,'\');
        anName=anFolders(anNum).name;
        thedate=anFolders(anNum).folder((thebreaks(2)+1):end);
        roiFile=[anName '_' thedate '_ROIs.mat'];
        cd(folderROIs);
        %% Load ROI file
        load(roiFile);
        roiNames=fieldnames(roiStruct);
        cd([anFolders(anNum).folder '\' anFolders(anNum).name]);
        VideoFiles=dir('*.avi');
       %% Process Video Files 
       allFileStart=tic;
        for filNum=1:size(VideoFiles,1)
            vidObj=VideoReader(VideoFiles(filNum).name); %Create an object containing video data and metadata
            trackCam.Params.frameWidth(filNum)=vidObj.Width;
            trackCam.Params.frameHeight(filNum)=vidObj.Height;
            trackCam.Params.frameCount(filNum)=vidObj.NumFrames;
            trackCam.Params.frameRate(filNum)=vidObj.FrameRate;
            trackCam.Params.timeStamp(filNum)=vidObj.CurrentTime;
            if filNum==1
                for roiNum=1:size(roiNames,1)
                    theDelta.(roiNames{roiNum})((1:trackCam.Params.frameCount(filNum)),(1:size(VideoFiles,1)))=NaN;
                    roiMat.(roiNames{roiNum})=repmat(uint8(roiStruct.(roiNames{roiNum}).roiMask),1,1,trackCam.Params.frameCount(filNum));
                    [z,p,k]=butter(3,5/(0.5*round(trackCam.Params.frameRate(filNum))),'low');
                    [sos_whisk,g_whisk]=zp2sos(z,p,k);
                end
            end

            %% Whole Video processing
            wholeFile=tic;
            videoLoad=tic;
            tempVid=read(vidObj,[1 Inf]);
            clear vidObj
            theVid=squeeze(tempVid(:,:,1,:));%Only use the first of the RGB channels and transform to HxWxT matrix
            clear tempVid
            video_load_time=toc(videoLoad)
            for roiNum=1:size(roiNames,1) 
                maskVid=theVid.*roiMat.(roiNames{roiNum})(:,:,(1:size(theVid,3)));
                diffVid=sum(sum(abs(diff(maskVid,1,3)),1),2);
                theDelta.(roiNames{roiNum})(1,filNum)=diffVid(1);
                theDelta.(roiNames{roiNum})((2:length(diffVid)+1),filNum)=diffVid;
                clear maskVid diffVid
            end
            Process_whole_trial_time=toc(wholeFile)
            %% Framewise processing IGNORE
%             for frameNum=2:trackCam.Params.frameCount(filNum)
%                 tic
%                 readTime=tic;
%                 leadFrame=read(vidObj,(frameNum-1));
%                 followFrame=read(vidObj,frameNum);
%                 readDur(frameNum)=toc(readTime);
%                 for roiNum=1:size(roiNames,1)  
%                         leadFramebin=double(leadFrame(:,:,1)).*double(roiStruct.(roiNames{roiNum}).roiMask);
%                         followFramebin=double(followFrame(:,:,1)).*double(roiStruct.(roiNames{roiNum}).roiMask);
%                     if frameNum==2
%                         theDelta.(roiNames{roiNum})((frameNum-1),filNum)=sum(abs(followFramebin-leadFramebin),'all');%duplicate first difference to prevent temporal shift for each file
%                     end
%                     theDelta.(roiNames{roiNum})((frameNum),filNum)=sum(abs(followFramebin-leadFramebin),'all');
%                 end
%                 frameTime(frameNum)=toc;
%             end
%             average_frame_read_time=mean(readDur)
%             average_frame_process_time=mean(frameTime)
%             total_file_process_time=sum(frameTime)
clear theVid
fprintf(['Finished tracking ' num2str(filNum) ' of ' num2str(size(VideoFiles,1)) ' video files.\n'])
        end
        total_process_time=toc(allFileStart)/60
        for roiNum=1:size(roiNames,1)
            catDelta=reshape(theDelta.(roiNames{roiNum}),1,numel(theDelta.(roiNames{roiNum}))); %concatenate video files in to single vector
            lastFrame=find(isnan(catDelta),1,'first')-1;
            smoothDelta=filtfilt(sos_whisk,g_whisk,double(catDelta(1:lastFrame))); %low-pass filter behavior data below 5Hz
            trackCam.Tracking.(roiNames{roiNum}).rawTracking=catDelta;
            trackCam.Tracking.(roiNames{roiNum}).smoothTracking=smoothDelta;
            figure;plot((1:length(smoothDelta))/trackCam.Params.frameRate(filNum),smoothDelta,'k','LineWidth',1);
            xlabel('Time (sec)'); ylabel('a.u.'); title([roiNames{roiNum} ' movement detection']);
        end
        clear theDelta
        cd(folderTracking);
        save([anName '_' thedate '_behaviorTracking.mat'],'trackCam','-v7.3'); %save movement tracking to be binarized later.
    end
end
end
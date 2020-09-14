function [binWhisk]=binarize_Whisking(VideoFolder,binThresh,varargin)
% This function loads in .avi files collected using BFLY usb cameras,
% allows the user to define an ROI around the whisker area and categorize
% whisking behavior in to a binary Whisking (1) or not Whisking(0) vector
% for analysis of imaging and electrophysiology data.

cd(VideoFolder);
VideoFiles=dir('*.avi');

for filNum=1:size(VideoFiles,1)
    VidObj=VideoReader(VideoFiles(filNum).name); %Create an object containing video data and metadata
    frameWidth=VidObj.Width;
    frameHeight=VidObj.Height;
    frameCount=VidObj.NumFrames;
    frameRate=round(VidObj.FrameRate);
    [z,p,k]=butter(3,5/(0.5*frameRate),'low');
    [sos_whisk,g_whisk]=zp2sos(z,p,k);
    if filNum==1
        FirstFrame=read(VidObj,1);
        whiskerFig=figure; whiskerAxes=axes(whiskerFig);
        whiskerImage=imshow(FirstFrame,'Parent',whiskerAxes);
        title(whiskerAxes,'Draw a rectangle around whiskers');
        whiskerROI=drawrectangle;
        confirm=input('Is ROI correct? (y/n)','s');
        if strcmpi(confirm,'y')
        theMask=createMask(whiskerROI);
        else
         whiskerROI=drawrectangle;
         theMask=createMask(whiskerROI);
        end
    end
    for frameNum=2:frameCount
        if frameNum==2
        leadFrame=read(VidObj,(frameNum-1));
        leadFramebin=double(leadFrame(:,:,1)).*theMask;
        else
        leadFramebin=followFramebin;
        end
        followFrame=read(VidObj,frameNum);
        followFramebin=double(followFrame(:,:,1)).*theMask;
        theDelta((frameNum-1),filNum)=sum(abs(followFramebin-leadFramebin),'all');
    end
end
smoothDelta=reshape(filtfilt(sos_whisk,g_whisk,theDelta),1,numel(theDelta));
figure;plot((1:length(smoothDelta))/frameRate,smoothDelta,'k','LineWidth',1);
xlabel('Time (sec)'); ylabel('a.u.'); title('Whisker movement detection');
if exist('binthresh','var')==0
    binThresh=input('Define threshold for whisker binarization');
end
binWhisk=smoothDelta>=binThresh;
end
function defineBehaviorROIs(VideoDirectory,roiDir)
%% READ ME
% Inputs:
%VideoDirectory: File directory containing behavior camera files to get ROI
% Directory format  established as 'DriveLetter:\Kyle\date\animal\files.avi'
%
%roiDir: File directory to save session specific ROI structures.


cd(VideoDirectory);
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
        cd([anFolders(anNum).folder '\' anFolders(anNum).name]);
        thebreaks=strfind(anFolders(anNum).folder,'\');
        anName=anFolders(anNum).name;
        thedate=anFolders(anNum).folder((thebreaks(2)+1):end);
        saveName=[anName '_' thedate '_ROIs.mat'];
        cd(roiDir);
        if ~exist(saveName, 'file')
            newROI='y';
        else
            newROI=input('ROI exists. Re-draw ROI? (y/n)','s');
        end
        if strcmpi(newROI,'y')
            cd([anFolders(anNum).folder '\' anFolders(anNum).name]);
            theFiles=dir('*.avi');
            vidObj=VideoReader(theFiles(1).name);
            FirstFrame=read(vidObj,1);
            confirm='n';
            roiNum=input('How many ROI do you want to define?');
            for roiInd=1:roiNum
                while(strcmpi(confirm,'n'))
                    whiskerFig=figure; whiskerAxes=axes(whiskerFig);
                    whiskerImage=imshow(FirstFrame,'Parent',whiskerAxes);
                    title(whiskerAxes,'Draw a rectangle around area you want to track');
                    theROI=drawrectangle;
                    confirm=input('Is ROI correct? (y/n)','s');
                    if strcmpi(confirm,'n')
                        close(whiskerFig);
                    end
                end
                roiName=input('What is the name of the ROI?', 's');
                roiMask=createMask(theROI);
                roiStruct.(roiName).roiMask=roiMask;
                roiStruct.(roiName).roiObj=theROI;
                close(whiskerFig);
                clear whiskerFig roiName roiMask;
                confirm='n';
            end
            cd(roiDir);
            saveName=[anName '_' thedate '_ROIs.mat'];
            save(saveName,'roiStruct');
        end
    end
end
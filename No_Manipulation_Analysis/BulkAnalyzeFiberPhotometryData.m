function BulkAnalyzeFiberPhotometryData
%READ ME
%This function will return a completed dataset from raw .csv files
%collected using Doric fiberphotometry data acquisition systems
%
%%Filenaming Format:
%%'ProjectIdentifier_AnimalNumber_Genotype_Fluorphore.csv'
%
%%Data Storage Format:
%%'DriveName\DoricData\AnimalGenotype\ProbeDepth\Filename.csv'%

%Embeded Functions
%Calibrate_Correction
%ExtractFPdata_003
%FiberGroupAnimal
%AssemblePopulationData_002
%CreatePlotStruct
%MakeDataTable
%GLME_FiberData_Fixed
%fiberphotometry_figures_004

Home='H:\QZ_FBR_FC_Data';
GFPDir='H:\DoricData\CAG_eGFP';
AnalogChannelNames={'RotaryEncoder'};
OpticalChannelNames={'Control','GCaMP6s','BloodVolume'};
popDir='H:\QZ_FBR_FC_Data\PopulationData';
folderTracking='H:\BehaviorCameraTracks';
%% Get correction for Hemodynamic attenuation of GFP signal
cd(GFPDir);
% if exist('Hemodynamic_Correction_Constant.mat','file')
%     load('Hemodynamic_Correction_Constant.mat');
% else
depthFolders=dir;
depthFolders(~[depthFolders.isdir])=[];
tf=ismember({depthFolders.name},{'.','..'});
depthFolders(tf)=[];
for foldNum=1:size(depthFolders)
    cd([depthFolders(foldNum).folder '\' depthFolders(foldNum).name]);
    TheFiles=dir(fullfile([depthFolders(foldNum).folder '\' depthFolders(foldNum).name],'**','*.csv'));
    for fileNum=1:size(TheFiles,1)
        if ~strcmpi(TheFiles(fileNum).folder,'E:\DoricData\CAG_eGFP\1500_um\082119\CE_FBR002')
            cd(TheFiles(fileNum).folder);
            filename=TheFiles(fileNum).name;
            [coeffVals,theEqn,goodness,stats,correctionConstant_xcorr,initCorrCoeffs,adjCorrCoeffs]=Calibrate_Correction(filename);
            FitData(fileNum).coeffVals=coeffVals;
            FitData(fileNum).fitEqn=theEqn;
            FitData(fileNum).GoodnessofFit=goodness;
            FitData(fileNum).Stats=stats;
            FitData(fileNum).minXcorr.correctionConst_xcorr=correctionConstant_xcorr;
            FitData(fileNum).minXcorr.initCorrCoeffs=initCorrCoeffs;
            FitData(fileNum).minXcorr.adjCorrCoeffs=adjCorrCoeffs;
            Slope(fileNum)=FitData(fileNum).coeffVals(1);
            Slope_xcorr(fileNum)=correctionConstant_xcorr;
        end
        if fileNum==size(TheFiles,1)
            Slope(Slope>=0)=NaN;
            CorrectionConst=nanmean(Slope);
            xcorrConst=mean(Slope_xcorr);
        end
    end
end
% cd(GFPDir);
% save('Hemodynamic_Correction_Constant','CorrectionConst');
% end
cd(Home);
subfolders=dir;
subfolders(~[subfolders.isdir])=[];
tf=ismember({subfolders.name},{'.','..','Ignore','PopulationData'});
subfolders(tf)=[];
%% Process raw .csv files
for folderNum=1:size(subfolders,1)
    cd([subfolders(folderNum).folder '\' subfolders(folderNum).name]);
    FileList=dir(fullfile([subfolders(folderNum).folder '\' subfolders(folderNum).name],'**','*.csv'));
    for filNum=1:size(FileList,1)
        cd(FileList(filNum).folder);
        filename=FileList(filNum).name;
        theBreaks=strfind(filename,'_');
        anName=filename(1:(theBreaks(2)-1));
        theSlash=strfind(FileList(filNum).folder,'\');
        filDate=FileList(filNum).folder((theSlash(end)+1):end);
        trackFile=[anName '_' filDate '_behaviorTracking.mat'];
        cd(folderTracking);
        if exist(trackFile,'file')
           trackData=load(trackFile);
        else
            fprintf('This file does not have behavior camera data\n')
            pause
        end
        cd(FileList(filNum).folder);
        correctFlag='n'; %If GCaMP signals need to be corrected for hemodynamic attenuation =='y' 
        ExtractFPdata_FC(filename,OpticalChannelNames,AnalogChannelNames,trackData.trackCam,correctFlag);
    end
end
close all;
% Group multiple imaging sessions for each unique animal
for folderNum=1:size(subfolders,1)
    cd([subfolders(folderNum).folder '\' subfolders(folderNum).name]);
    FileList=dir(fullfile([subfolders(folderNum).folder '\' subfolders(folderNum).name],'**','QZ_FBR*.mat'));
    for filNum=1:size(FileList,1)
        thebreaks=strfind(FileList(filNum).name,'_');
        Animal_Name{filNum}=FileList(filNum).name(1:(thebreaks(2)-1));
    end
    Ind_Animals=unique(Animal_Name);
    for anNum=1:size(Ind_Animals,2)
        anFind=strcmpi(Animal_Name,Ind_Animals{anNum});
        theFiles=FileList(anFind);
        FiberGroupAnimal(theFiles,popDir);
    end   
    clear Animal_Name
end

%% Average animals in each subgroup
cd(popDir);
FileList=dir(fullfile(popDir,'**','QZ_FBR*.mat'));
for filnum=1:size(FileList,1)
    name_break=strfind(FileList(filnum).name,'_');
    AnimalGroup{filnum}=FileList(filnum).name((name_break(2)+1):(name_break(4)-1));
end
Promoters=unique(AnimalGroup);
for proNum=1:size(Promoters,2)
    filFind=strcmpi(AnimalGroup,Promoters{proNum});
    TempAnimals=FileList(filFind);
    for anNum=1:size(TempAnimals,1)
        thebreaks=strfind(TempAnimals(anNum).name,'_');
        fiberDepth{anNum}=TempAnimals(anNum).name((thebreaks(4)+1):(thebreaks(5)-1));%Changed end to 5 from 6 for QZ project
    end
    Depths=unique(fiberDepth);
    if proNum==1
        GroupedData=[];
    end
   
    theAnimals=TempAnimals;
    [GroupedData]=AssemblePopulationData_002(theAnimals,popDir,Promoters{proNum},Depths{1},GroupedData);
    clear fiberDepth
end
clear AnimalGroup

%% Save Data
saveDate=datetime('now');
theDate=datestr(saveDate);
theDate=strrep(theDate,'-','');
theDate=strrep(theDate,' ','_');
theDate=strrep(theDate,':','');
cd(popDir);
if exist('Fitdata','var')
save(['GroupedData_' theDate '_FiberPhotometry'],'GroupedData','FitData','-v7.3');
else
save(['GroupedData_' theDate '_FiberPhotometry'],'GroupedData','-v7.3');  
end

%% Create Data structure for figures
[FigureData]=CreatePlotStruct(GroupedData); % Need to adapt for QZ Project

%% Create/Apply Generalized Linear Model between datasets
[Locomotion_Data_Table]=MakeDataTable_002(GroupedData,FigureData);
[GLME_Output]=GLME_FiberData_Fixed(Locomotion_Data_Table);
FigureData.Stats.StatsData=Locomotion_Data_Table;
FigureData.Stats.GLME_Fits=GLME_Output;

saveDate=datetime('now');
theDate=datestr(saveDate);
theDate=strrep(theDate,'-','');
theDate=strrep(theDate,' ','_');
theDate=strrep(theDate,':','');
save(['FigureData_' theDate],'FigureData','-v7.3');
%% Visualize Fiber photometry data
fiberphotometry_figures_004(FigureData);
end
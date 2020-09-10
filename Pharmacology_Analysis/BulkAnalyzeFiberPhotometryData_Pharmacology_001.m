function BulkAnalyzeFiberPhotometryData_Pharmacology_001

Home='H:\DREADDs_FiberPhotometry';
GFPDir='H:\DoricData\CAG_eGFP';
AnalogChannelNames={'RotaryEncoder'};
OpticalChannelNames={'Control','GCaMP6s','TRITC'};
popDir='H:\DREADDs_FiberPhotometry\PopulationData';

%% Get correction for Hemodynamic attenuation of GFP signal
cd(GFPDir);
depthFolders=dir;
depthFolders(~[depthFolders.isdir])=[];
tf=ismember({depthFolders.name},{'.','..'});
depthFolders(tf)=[];
for foldNum=1:size(depthFolders)
    cd([depthFolders(foldNum).folder '\' depthFolders(foldNum).name]);
    TheFiles=dir(fullfile([depthFolders(foldNum).folder '\' depthFolders(foldNum).name],'**','*.csv'));
    for fileNum=1:size(TheFiles,1)
        if ~strcmpi(TheFiles(fileNum).folder,'H:\DoricData\CAG_eGFP\1500_um\082119\CE_FBR002')
            cd(TheFiles(fileNum).folder);
            filename=TheFiles(fileNum).name;
            [coeffVals,theEqn,goodness,stats]=Calibrate_Correction_Pharmacology(filename);
            FitData(fileNum).coeffVals=coeffVals;
            FitData(fileNum).fitEqn=theEqn;
            FitData(fileNum).GoodnessofFit=goodness;
            FitData(fileNum).Stats=stats;
            Slope(fileNum)=FitData(fileNum).coeffVals(1);
        end
        if fileNum==size(TheFiles,1)
            Slope(Slope>=0)=NaN;
            CorrectionConst=nanmean(Slope);
        end
    end
end
close all

%% Process raw .csv files
cd(Home);
subfolders=dir;
subfolders(~[subfolders.isdir])=[];
tf=ismember({subfolders.name},{'.','..','Ignore','PopulationData'});
subfolders(tf)=[];

for folderNum=1:size(subfolders,1)
    cd([subfolders(folderNum).folder '\' subfolders(folderNum).name]);
    FileList=dir(fullfile([subfolders(folderNum).folder '\' subfolders(folderNum).name],'**','*.csv'));
    
    for filNum=1:size(FileList,1)
        animalLabel{filNum}=FileList(filNum).name(1:9);
    end
    animalName=unique(animalLabel);
    
    for anNum=1:size(animalName,2)
        filename=[];
        fileCnt=1;
        for filNum=1:size(FileList,1)
            anFiles(filNum)=contains(FileList(filNum).name,animalName{anNum},'IgnoreCase',true);
            if anFiles(filNum)==1
                folderLoc=FileList(filNum).folder;
                filename{fileCnt}=FileList(filNum).name;
                writeTime(fileCnt)=datetime(FileList(filNum).date);
                fileCnt=fileCnt+1;
            end
        end
        
        cd(folderLoc);
        ExtractFPdata_Pharmacology_001(filename,OpticalChannelNames,AnalogChannelNames,CorrectionConst,writeTime);
        clear writeTime
    end
end
close all;

%% Group multiple imaging sessions for each unique animal
cd(subfolders(1).folder);
FileList=dir(fullfile(subfolders(1).folder,'**','CE_FBR*FiberPhotometry.mat'));
for filNum=1:size(FileList,1)
    thebreaks=strfind(FileList(filNum).name,'_');
    Animal_Name{filNum}=FileList(filNum).name(1:(thebreaks(2)-1));
end
Ind_Animals=unique(Animal_Name);
for anNum=1:size(Ind_Animals,2)
    anFind=strcmpi(Animal_Name,Ind_Animals{anNum});
    theFiles=FileList(anFind);
    for filNum=1:size(theFiles,1)
        theBreaks=strfind(theFiles(filNum).folder,'\');
        Treatment{filNum}=theFiles(filNum).folder((theBreaks(3)+1):(theBreaks(4)-1));
    end
    TreatmentType=unique(Treatment);
    for treatNum=1:size(TreatmentType,2)
        treatFind=strcmpi(Treatment,TreatmentType{treatNum});
        procFiles=theFiles(treatFind);
        FiberGroupAnimal_Pharmacology(procFiles,popDir);
    end
end
clear Animal_Name

%% Average animals in each subgroup
cd(popDir);
FileList=dir(fullfile(popDir,'**','CE_FBR*.mat'));
for filnum=1:size(FileList,1)
    name_break=strfind(FileList(filnum).name,'_');
    Injection{filnum}=FileList(filnum).name((name_break(4)+1):(name_break(5)-1));

end
injectionTypes=unique(Injection);
for injNum=1:size(injectionTypes,2)
    filFind=strcmpi(Injection,injectionTypes{injNum});
    TempAnimals=FileList(filFind);
    for anNum=1:size(TempAnimals,1)
        thebreaks=strfind(TempAnimals(anNum).name,'_');
        AnimalGroup{anNum}=TempAnimals(anNum).name((thebreaks(2)+1):(thebreaks(4)-1));
    end
    Genotypes=unique(AnimalGroup);
    for genNum=1:size(Genotypes,2)
        filFind=strcmpi(AnimalGroup,Genotypes{genNum});
        genotypeAnimals=TempAnimals(filFind);
        for anNum=1:size(genotypeAnimals,1)
            thebreaks=strfind(genotypeAnimals(anNum).name,'_');
            fiberDepth{anNum}=genotypeAnimals(anNum).name((thebreaks(5)+1):(thebreaks(6)-1));
        end
    Depths=unique(fiberDepth);
    for depthNum=1:size(Depths,2)
    if injNum==1
        GroupedData=[];
    end
    filFind=strcmpi(fiberDepth,Depths{depthNum});
    theAnimals=genotypeAnimals(filFind);
    [GroupedData]=AssemblePopulationData_Pharmacology(theAnimals,popDir,injectionTypes{injNum},Genotypes{genNum},Depths{depthNum},GroupedData);
    end
    clear fiberDepth
    clear genotypeAnimals
    clear theAnimals
    end
    clear AnimalGroup
    clear TempAnimals
end
clear Injection

%% Save Data
saveDate=datetime('now');
theDate=datestr(saveDate);
theDate=strrep(theDate,'-','');
theDate=strrep(theDate,' ','_');
theDate=strrep(theDate,':','');
cd(popDir);
save(['GroupedData_' theDate '_FiberPhotometry'],'GroupedData','-v7.3');%,'FitData'

%% Create Data structure for figures
[FigureData]=CreatePlotStruct(GroupedData);

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
    
    
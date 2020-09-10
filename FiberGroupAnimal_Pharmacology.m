function FiberGroupAnimal_Pharmacology(Filedir,popDir)
%Filedir is a structure output from function 'dir' with file names and
%folder locations.

for fileNum=1:size(Filedir,1)
    cd(Filedir(fileNum).folder);
    thebreaks=strfind(Filedir(fileNum).folder,'\');
    load(Filedir(fileNum).name);
    animalName=ChunkData.Params.animalname;
    proName=ChunkData.Params.AnimalType;
    fiberDepth=ChunkData.Params.fiber_depth;
    injectionType=Filedir(fileNum).folder((thebreaks(3)+1):(thebreaks(4)-1));
    datafields=fieldnames(ChunkData);
    if fileNum==1
        AnimalData.Params=ChunkData.Params;
        AnimalAvg.Params=ChunkData.Params;
    end
    keepfields=ismember(datafields,{'Params','WheelData','LocomotionEvokedData','TrialData'});%'RawFiberData',
    datafields(keepfields)=[];
    for fieldNum=1:size(datafields,1)
        if strcmpi(datafields{fieldNum},'RawFiberData')
            AnimalData.(datafields{fieldNum}).RawSignalAvg(fileNum,:)=mean(ChunkData.RawFiberData.RawFiberData,1);
            AnimalData.(datafields{fieldNum}).DetrendSignalSTD(fileNum,:)=std(ChunkData.RawFiberData.RawFiberData,0,1);
            AnimalData.(datafields{fieldNum}).LowPassSignalSTD(fileNum,:)=std(ChunkData.RawFiberData.LowPassData,0,1);
            AnimalData.(datafields{fieldNum}).RescaleAvg(fileNum,:)=mean(ChunkData.RawFiberData.RescaledData,1);
            AnimalData.(datafields{fieldNum}).RescaleStd(fileNum,:)=std(ChunkData.RawFiberData.RescaledData,0,1);
            AnimalData.(datafields{fieldNum}).Corrected465Avg(fileNum,:)=mean(ChunkData.RawFiberData.Corrected465);
            AnimalData.(datafields{fieldNum}).Corrected465STD(fileNum,:)=std(ChunkData.RawFiberData.Corrected465,0,1);
            if fileNum==size(Filedir,1)
                AnimalAvg.(datafields{fieldNum}).RawSignalAvg=mean(AnimalData.(datafields{fieldNum}).RawSignalAvg,1);
                AnimalAvg.(datafields{fieldNum}).DetrendSignalSTD=mean(AnimalData.(datafields{fieldNum}).DetrendSignalSTD,1);
                AnimalAvg.(datafields{fieldNum}).LowPassSignalSTD=mean(AnimalData.(datafields{fieldNum}).LowPassSignalSTD,1);
                AnimalAvg.(datafields{fieldNum}).RescaleAvg=mean(AnimalData.(datafields{fieldNum}).RescaleAvg,1);
                AnimalAvg.(datafields{fieldNum}).RescaleStd=mean(AnimalData.(datafields{fieldNum}).RescaleStd,1);
                AnimalAvg.(datafields{fieldNum}).Corrected465Avg=mean(AnimalData.(datafields{fieldNum}).Corrected465Avg,1);
                AnimalAvg.(datafields{fieldNum}).Corrected465STD=mean(AnimalData.(datafields{fieldNum}).Corrected465STD,1);
            end
        else
            subfields=fieldnames(ChunkData.(datafields{fieldNum}));
            for subNum=1:size(subfields,1)
                if isstruct(ChunkData.(datafields{fieldNum}).(subfields{subNum}))
                    finalNames=fieldnames(ChunkData.(datafields{fieldNum}).(subfields{subNum}));
                    for finNum=1:size(finalNames,1)
                        if isstruct(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}))
                            lastNames=fieldnames(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}));
                            for lastNum=1:size(lastNames,1)
                                if ~isempty(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum}))
                                    if fileNum~=1
                                        if length(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum})(:,:,(fileNum-1)))==1
                                            AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum})=[];
                                            AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum})=NaN(size(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum}),1),...
                                                size(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum}),2),(fileNum-1));
                                        end
                                        AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum})(:,:,fileNum)=ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum});
                                    else
                                        AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum})(:,:,fileNum)=ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum});
                                    end
                                else
                                    if fileNum==1
                                        AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum})(:,:,fileNum)=nan;
                                    else
                                        AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum})(:,:,fileNum)=NaN(size(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum})(:,:,(fileNum-1)),1),...
                                            size(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum})(:,:,(fileNum-1)),2));
                                    end
                                    
                                end
                                
                                if fileNum==size(Filedir,1)
                                    AnimalAvg.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum})=nanmean(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{lastNum}),3);
                                end
                            end
                        else
                            if ~isempty(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}))
                                if fileNum~=1
                                    if length(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum})(:,:,(fileNum-1)))==1
                                        AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum})=[];
                                        AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum})=NaN(size(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}),1),size(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}),2),(fileNum-1));
                                    end
                                    AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum})(:,:,fileNum)=ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum});
                                else
                                    AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum})(:,:,fileNum)=ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum});
                                end
                            else
                                if fileNum==1
                                    AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum})(:,:,fileNum)=nan;
                                else
                                    AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum})(:,:,fileNum)=NaN(size(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum})(:,:,(fileNum-1)),1),...
                                        size(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum})(:,:,(fileNum-1)),2));
                                end
                                
                            end
                            
                            if fileNum==size(Filedir,1)
                                AnimalAvg.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum})=nanmean(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}),3);
                            end
                        end
                    end
                else
                    AnimalData.(datafields{fieldNum}).(subfields{subNum})(:,:,fileNum)=ChunkData.(datafields{fieldNum}).(subfields{subNum});
                    if fileNum==size(Filedir,1)
                        AnimalAvg.(datafields{fieldNum}).(subfields{subNum})=nanmean(AnimalData.(datafields{fieldNum}).(subfields{subNum}),3);
                    end
                end
            end
        end
    end
end
cd(popDir);
save([animalName '_' proName '_' injectionType '_' fiberDepth '_AnimalAvgData'],'AnimalData','AnimalAvg','-v7.3');
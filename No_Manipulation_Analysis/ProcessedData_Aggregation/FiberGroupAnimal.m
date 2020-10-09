function FiberGroupAnimal(Filedir,popDir)
%Filedir is a structure output from function 'dir' with file names and
%folder locations.

for fileNum=1:size(Filedir,1)
    cd(Filedir(fileNum).folder);
    thebreaks=strfind(Filedir(fileNum).folder,'\');
    load(Filedir(fileNum).name);
    animalName=ChunkData.Params.animalname;
    proName=Filedir(fileNum).folder((thebreaks(2)+1):(thebreaks(3)-1));
    fiberDepth=ChunkData.Params.fiber_depth;
    datafields=fieldnames(ChunkData);
    if fileNum==1
        AnimalData.Params=ChunkData.Params;
        AnimalAvg.Params=ChunkData.Params;
    end
    keepfields=ismember(datafields,{'Params','WheelData','LocomotionEvokedData'});%'RawFiberData',
    datafields(keepfields)=[];
    for fieldNum=1:size(datafields,1)
        if ~strcmpi(datafields{fieldNum},'WhiskersEvokedData') && ~strcmpi(datafields{fieldNum},'FaceEvokedData')
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
                            if ~isstruct(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}))
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
                            else
                                lastNames=fieldnames(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}));
                                for nameNum=1:size(lastNames,1)
                                    if ~isstruct(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}))
                                        if ~isempty(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}))
                                            if fileNum~=1
                                                if ~isfield(AnimalData.(datafields{fieldNum}).(subfields{subNum}),finalNames{finNum})
                                                    AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum})=[];
                                                end
                                                if ~isfield(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}),lastNames{nameNum})
                                                    AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum})=[];
                                                end
                                                structsize=size(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}),3);
                                                if length(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum})(:,:,structsize))==1
                                                    AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum})=[];
                                                    AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum})=NaN(size(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}),1),size(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}),2),(fileNum-1));
                                                end
                                                AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum})(:,:,fileNum)=ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum});
                                            else
                                                AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum})(:,:,fileNum)=ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum});
                                            end
                                        else
                                            if fileNum==1
                                                AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum})(:,:,fileNum)=nan;
                                            else
                                                if ~isfield(AnimalData.(datafields{fieldNum}).(subfields{subNum}),finalNames{finNum})
                                                    AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum})=[];
                                                end
                                                if ~isfield(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}),lastNames{nameNum})
                                                    AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum})=[];
                                                end
                                                structsize=size(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}),3);
                                                AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum})(:,:,fileNum)=NaN(size(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum})(:,:,structsize),1),...
                                                    size(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum})(:,:,structsize),2));
                                            end
                                            
                                        end
                                        
                                        if fileNum==size(Filedir,1)
                                            AnimalAvg.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum})=nanmean(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}),3);
                                        end
                                    else
                                        termName=fieldnames(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}));
                                        for termNum=1:size(termName,1)
                                            if ~isempty(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum}))
                                                if fileNum~=1
                                                    if ~isfield(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}),lastNames{nameNum})
                                                        AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum})=[];
                                                    end
                                                    if ~isfield(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}),termName{termNum})
                                                        AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum})=[];
                                                    end
                                                    if length(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum})(:,:,(fileNum-1)))==1
                                                        AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum})=[];
                                                        AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum})=NaN(size(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum}),1),size(ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum}),2),(fileNum-1));
                                                    end
                                                    AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum})(:,:,fileNum)=ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum});
                                                else
                                                    AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum})(:,:,fileNum)=ChunkData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum});
                                                end
                                            else
                                                if fileNum==1
                                                    AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum})(:,:,fileNum)=nan;
                                                else
                                                    AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum})(:,:,fileNum)=NaN(size(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum})(:,:,(fileNum-1)),1),...
                                                        size(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum})(:,:,(fileNum-1)),2));
                                                end
                                                
                                            end
                                            
                                            if fileNum==size(Filedir,1)
                                                AnimalAvg.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum})=nanmean(AnimalData.(datafields{fieldNum}).(subfields{subNum}).(finalNames{finNum}).(lastNames{nameNum}).(termName{termNum}),3);
                                            end
                                        end  
                                    end
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
end
cd(popDir);
save([animalName '_' proName '_' fiberDepth '_AnimalAvgData'],'AnimalData','AnimalAvg','-v7.3');
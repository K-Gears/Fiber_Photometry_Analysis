function [GroupedData]=AssemblePopulationData_002(AnimalFiles,popDir,proName,DepthName,GroupedData,varargin)
if ~exist('GroupedData','var')
GroupedFile=dir(fullfile(popDir,'**','Grouped*.mat')); %check for existence of grouped data file
if ~isempty(GroupedFile)
    cd(GroupedFile.folder);
    load(GroupedFile.name);
else
    GroupedData=[];
end
end
depthBreak=strfind(DepthName,'_');
depthVal=str2double(DepthName(1:(depthBreak(1)-1)));
if depthVal>=1200
    DepthName='Deep';
elseif depthVal>=750
    DepthName='Mid';
else
    DepthName='Shallow';
end

for anNum=1:size(AnimalFiles,1)
    load(AnimalFiles(anNum).name,'AnimalAvg');
    nameBreak=strfind(AnimalFiles(anNum).name,'_');
    animalName=AnimalFiles(anNum).name(1:(nameBreak(2)-1));
    subfields=fieldnames(AnimalAvg);
    if anNum==1
        GroupedData.(proName).Params=AnimalAvg.Params;
    end
    GroupedData.(proName).Params.animalNames{anNum}=animalName;
    keepfields=~ismember(subfields,{'Params','Averages'});
    subfields=subfields(keepfields);
    for subNum=1:size(subfields,1)
        datafields=fieldnames(AnimalAvg.(subfields{subNum}));
        for dataNum=1:size(datafields,1)
            if isstruct(AnimalAvg.(subfields{subNum}).(datafields{dataNum}))
                finalfields=fieldnames(AnimalAvg.(subfields{subNum}).(datafields{dataNum}));
                for fieldNum=1:size(finalfields,1)
                    if anNum>1
                        if isnan(GroupedData.(proName).(subfields{subNum}).(datafields{dataNum}).(finalfields{fieldNum})(:,:,(anNum-1)))
                            GroupedData.(proName).(subfields{subNum}).(datafields{dataNum}).(finalfields{fieldNum})(:,:,(anNum-1))=[];
                            GroupedData.(proName).(subfields{subNum}).(datafields{dataNum}).(finalfields{fieldNum})((1:size(AnimalAvg.(subfields{subNum}).(datafields{dataNum}).(finalfields{fieldNum}),1)),(1:size(AnimalAvg.(subfields{subNum}).(datafields{dataNum}).(finalfields{fieldNum}),2)),(anNum-1))...
                                =NaN;
                        end
                    end
                    GroupedData.(proName).(subfields{subNum}).(datafields{dataNum}).(finalfields{fieldNum})(:,:,anNum)=AnimalAvg.(subfields{subNum}).(datafields{dataNum}).(finalfields{fieldNum});
                    if anNum==size(AnimalFiles,1)
                        if strcmpi(subfields{subNum},'RawFiberData')
                            GroupedData.(proName).Averages.(subfields{subNum}).(datafields{dataNum}).(finalfields{fieldNum})=median(GroupedData.(proName).(subfields{subNum}).(datafields{dataNum}).(finalfields{fieldNum}),3);
                        else
                        GroupedData.(proName).Averages.(subfields{subNum}).(datafields{dataNum}).(finalfields{fieldNum})=mean(GroupedData.(proName).(subfields{subNum}).(datafields{dataNum}).(finalfields{fieldNum}),3);
                        end
                    end
                end
            else
                GroupedData.(proName).(subfields{subNum}).(datafields{dataNum})(:,:,anNum)=AnimalAvg.(subfields{subNum}).(datafields{dataNum});
                if anNum==size(AnimalFiles,1)
                    if strcmpi(subfields{subNum},'RawFiberData')
                        GroupedData.(proName).Averages.(subfields{subNum}).(datafields{dataNum})=median(GroupedData.(proName).(subfields{subNum}).(datafields{dataNum}),3);
                    else
                        GroupedData.(proName).Averages.(subfields{subNum}).(datafields{dataNum})=mean(GroupedData.(proName).(subfields{subNum}).(datafields{dataNum}),3);
                    end
                end
            end
        end
    end
end
end
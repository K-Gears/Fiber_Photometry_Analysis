function [GroupedFPData]=GroupFiberPhotometry
GroupedFPData.LocomotionData.AveragedData=[];
GroupedFPData.LocomotionData.StdData=[];
for filenum=1:size(TheFiles,1)
    load(TheFiles(filenum).name);
    if filenum==1
        GroupedFPData.Params=ChunkData.Params;
    end
    GroupedFPData.LocomotionData.AveragedData=cat(3,GroupedFPData.LocomotionData.AveragedData,ChunkData.AveragedData.AvgResp);
    GroupedFPData.LocomotionData.StdData=cat(3,GroupedFPData.LocomotionData.StdData,ChunkData.AveragedData.StdResp);
    GroupedFPData.LocomotionData.EventCounts(filenum)=size(ChunkData.LocomotionData.OpticalData,3);
end
GroupedFPData.PopulationAverages.Locomotion=mean(GroupedFPData.LocomotionData.AveragedData,3);
animalnum=num2str(size(GroupedFPData.LocomotionData.AveragedData,3));
figure(99)
LocoTime=((1:length(GroupedFPData.PopulationAverages.Locomotion))-GroupedFPData.Params.StartPad)/GroupedFPData.Params.DataFs;
plot(LocoTime,GroupedFPData.PopulationAverages.Locomotion);
title(['Locomotion evoked change (n= ' animalnum ' )']);
xlabel('Time (s)');
ylabel('Z-score');
xlim([-3 5]);
lgndtxt={'405nm Ca2+ Independent','465nm Ca2+ dependent','560nm TRITC Blood Volume'};
legend(lgndtxt);

figure(100)
AnNum=1;
trialcount=num2str(GroupedFPData.LocomotionData.EventCounts(AnNum));
LocoTime=((1:length(GroupedFPData.PopulationAverages.Locomotion))-GroupedFPData.Params.StartPad)/GroupedFPData.Params.DataFs;
plot(LocoTime,GroupedFPData.LocomotionData.AveragedData(:,:,AnNum)); hold on;
for stanNum=1:size(GroupedFPData.LocomotionData.AveragedData,2)
    PosStan=GroupedFPData.LocomotionData.AveragedData(:,stanNum,AnNum)+GroupedFPData.LocomotionData.StdData(:,stanNum,AnNum);
    NegStan=GroupedFPData.LocomotionData.AveragedData(:,stanNum,AnNum)-GroupedFPData.LocomotionData.StdData(:,stanNum,AnNum);
    plot(LocoTime,PosStan,'--');plot(LocoTime,NegStan,'--');
end
title(['Locomotion evoked change (n= ' trialcount ' )']);
xlabel('Time (s)');
ylabel('Z-score');
xlim([-3 5]);
lgndtxt={'405nm Ca2+ Independent','465nm Ca2+ dependent','560nm TRITC Blood Volume'};
legend(lgndtxt);
save('FiberPhotometryPopulationData.mat','GroupedFPData','-nocompression','-v7.3');
    
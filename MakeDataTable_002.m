function [Locomotion_Data_Table]=MakeDataTable_002(GroupedData,FigureData)
AnimalNames=[];
AnimalGenotype=[];
AnimalProbeDepth=[];
CBVPeak5s=[];
GCaMPPeak5s=[];
CBVPeak10s=[];
GCaMPPeak10s=[];
CBVPeak15s=[];
GCaMPPeak15s=[];
CBVPeak30s=[];
GCaMPPeak30s=[];
XCorrCoeffPeak=[];
XCorrCoeffLag=[];
RawGCaMPAvg=[];
RawGCaMPStd=[];
RawCBVAvg=[];
RawCBVStd=[];

popDir='H:\DoricData\PopulationData';
cd(popDir);
Genotypes=fieldnames(GroupedData);

for folderNum=2:size(Genotypes)
    AnimalNames=[AnimalNames;GroupedData.(Genotypes{folderNum}).Params.animalNames'];
    for anNum=1:size(GroupedData.(Genotypes{folderNum}).Params.animalNames,2)
        SearchHandle=[GroupedData.(Genotypes{folderNum}).Params.animalNames{anNum}, '*.mat'];
        FileList=dir(fullfile(popDir,'**',SearchHandle));
        txtBreaks=strfind(FileList.name,'_');
        ProbeDepth(anNum)=str2double(FileList.name((txtBreaks(4)+1):(txtBreaks(5)-1)));
        TempType{anNum}=Genotypes{folderNum};
    end
    AnimalGenotype=[AnimalGenotype;TempType'];
    AnimalProbeDepth=[AnimalProbeDepth;ProbeDepth'];
    clear TempType ProbeDepth;
    
    TempData5s=squeeze(GroupedData.(Genotypes{folderNum}).AveragedData.five_second_events.PeakResp);
    CBVPeak5s=[CBVPeak5s;TempData5s(3,:)'];
    GCaMPPeak5s=[GCaMPPeak5s;TempData5s(2,:)'];
    
    TempData10s=squeeze(GroupedData.(Genotypes{folderNum}).AveragedData.ten_second_events.PeakResp);
    CBVPeak10s=[CBVPeak10s;TempData10s(3,:)'];
    GCaMPPeak10s=[GCaMPPeak10s;TempData10s(2,:)'];
    
    TempData15s=squeeze(GroupedData.(Genotypes{folderNum}).AveragedData.fifteen_second_events.PeakResp);
    CBVPeak15s=[CBVPeak15s;TempData15s(3,:)'];
    GCaMPPeak15s=[GCaMPPeak15s;TempData15s(2,:)'];
    
    TempData30s=squeeze(GroupedData.(Genotypes{folderNum}).AveragedData.thirty_second_events.PeakResp);
    CBVPeak30s=[CBVPeak30s;TempData30s(3,:)'];
    GCaMPPeak30s=[GCaMPPeak30s;TempData30s(2,:)'];
    
    XCorrCoeffPeak=[XCorrCoeffPeak;FigureData.Xcorr.CrossCorr.(Genotypes{folderNum}).PeakVal'];
    
    XCorrCoeffLag=[XCorrCoeffLag;FigureData.Xcorr.CrossCorr.(Genotypes{folderNum}).PeakTime'];
    
    RawGCaMPAvg=[RawGCaMPAvg;FigureData.Xcorr.SNR.(Genotypes{folderNum}).RawGCaMPAvg'];
    
    RawGCaMPStd=[RawGCaMPStd;log(FigureData.Xcorr.SNR.(Genotypes{folderNum}).RawGCaMPSTD')];
    
    RawCBVAvg=[RawCBVAvg;FigureData.Xcorr.SNR.(Genotypes{folderNum}).RawCBVAvg'];
    
    RawCBVStd=[RawCBVStd;log(FigureData.Xcorr.SNR.(Genotypes{folderNum}).RawCBVSTD')];
    
end

Locomotion_Data_Table=table(AnimalNames,AnimalGenotype,AnimalProbeDepth,GCaMPPeak5s,CBVPeak5s,GCaMPPeak10s,...
    CBVPeak10s,GCaMPPeak15s,CBVPeak15s,GCaMPPeak30s,CBVPeak30s,XCorrCoeffPeak,XCorrCoeffLag,RawGCaMPAvg,RawGCaMPStd,RawCBVAvg,RawCBVStd);

save('DataforGLME','Locomotion_Data_Table','-v7.3');
end
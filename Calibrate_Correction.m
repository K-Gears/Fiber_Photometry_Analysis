function [coeffVals,theEqn,goodness,stats]=Calibrate_Correction(filename)
%READ ME
%This function estimates the attenuation of eGFP signals due to increases
%in blood volume using data recorded from CAG-EGFP mice in L5 of S1 cortex.

close all
directory=cd;
folderBreaks=strfind(directory,'\');
ChunkData.Params.animalname=filename(1:9);
if numel(folderBreaks)==4
    ChunkData.Params.date=directory((folderBreaks(4)+1):end);
    ChunkData.Params.GFP_Type=directory((folderBreaks(2)+1):(folderBreaks(3)-1));
    ChunkData.Params.fiber_depth=directory((folderBreaks(3)+1):(folderBreaks(4)-1));
    
else
    ChunkData.Params.date=directory((folderBreaks(4)+1):(folderBreaks(5)-1));
    ChunkData.Params.GFP_Type=directory((folderBreaks(2)+1):(folderBreaks(3)-1));
    ChunkData.Params.fiber_depth=directory((folderBreaks(3)+1):(folderBreaks(4)-1));
end
ChunkData.Params.Acquisition_Fs=1.2e4;
ChunkData.Params.Decimation=10;
ChunkData.Params.DataFs=ChunkData.Params.Acquisition_Fs/ChunkData.Params.Decimation;
ChunkData.Params.DataAcquired=2;%time in hours of data acquired
ChunkData.Params.DataSeconds=ChunkData.Params.DataAcquired*(60^2);
ChunkData.Params.VelocityChannel=7;
ChunkData.Params.StartPad=5*ChunkData.Params.DataFs;%Time to collect in front of running start
ChunkData.Params.FollowPad=15*ChunkData.Params.DataFs;%Time to follow running start
ChunkData.Params.Fit_Freq=0.01;
ChunkData.Params.low_Freq=1;
ChunkData.Params.Final_Freq=[0.01 1];
params.Fs=ChunkData.Params.DataFs; %multitaper estimation parameters
params.tapers=[3 5];%multitaper estimation parameters

[z,p,k]=butter(3,ChunkData.Params.Fit_Freq/(0.5*ChunkData.Params.DataFs),'low'); %design lowpass filter for hemodynamic correction
[sos_Fit,g_Fit]=zp2sos(z,p,k);

[z,p,k]=butter(3,ChunkData.Params.low_Freq/(0.5*ChunkData.Params.DataFs),'low'); %Low pass for optical data to physiologically relevant range
[sos_Low,g_Low]=zp2sos(z,p,k);

[z,p,k]=butter(3,ChunkData.Params.Final_Freq/(0.5*ChunkData.Params.DataFs),'bandpass'); %Low pass filter for locomotion data
[sos_final,g_final]=zp2sos(z,p,k);

[z,p,k]=butter(3,10/(0.5*ChunkData.Params.DataFs),'low'); %Low pass filter for locomotion data
[sos_ball,g_ball]=zp2sos(z,p,k);

%% Read .CSV
fileInf=detectImportOptions(filename);
DataChannels=[2,3,5];%these are the three demodulated optical channels [2,3,5] and the wheel data [7]
RawData=csvread(filename,2,0);
WheelData=filtfilt(sos_ball,g_ball,detrend(RawData(:,ChunkData.Params.VelocityChannel))); %(5*ChunkData.Params.DataFs):(end-(5*ChunkData.Params.DataFs))
RawData=RawData(:,DataChannels); %(5*ChunkData.Params.DataFs):(end-(5*ChunkData.Params.DataFs))

%% Find Locomotion points to exclude from baseline calculations
[imp_bin]=velocity_binarize_fiberphotometry(WheelData,ChunkData.Params.DataFs,ChunkData.Params.DataFs,1e-3);
FuseGaps=15*ChunkData.Params.DataFs;
bin_Run=double(imp_bin);
RunInds=find(bin_Run==1);
IndGap=diff(RunInds);
for IndNum=1:length(IndGap)
    if IndGap(IndNum)>1
        if IndGap(IndNum)<=FuseGaps
            bin_Run(RunInds(IndNum):RunInds(IndNum+1))=1;
        else
            bin_Run(RunInds(IndNum):(RunInds(IndNum)+FuseGaps))=1;
        end
    end
end
ExcludeVals=find(bin_Run==1);

%% Remove photobleaching/metaloism baseline drift
Spacing=1:1:length(RawData(:,3));
FiltData=filtfilt(sos_Fit,g_Fit,RawData); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
%Correct TRITC blood volume
[fitVals]=fit(Spacing',FiltData(:,3),'exp2','Exclude',ExcludeVals);
coeffVals=coeffvalues(fitVals);
predictedCBV=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing))); 
CorrectedCBV=RawData(:,3)-predictedCBV';
FitStruct.CBV=fitVals;
% Correct Ca2+ dependent GCaMP
[fitVals]=fit(Spacing',FiltData(:,2),'exp2','Exclude',ExcludeVals);
coeffVals=coeffvalues(fitVals);
predictedGCaMP=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing))); 
CorrectedGCaMP=RawData(:,2)-predictedGCaMP';
FitStruct.GCaMP=fitVals;
% Correct Ca2+ independent GCaMP
[fitVals]=fit(Spacing',FiltData(:,1),'exp2','Exclude',ExcludeVals);
coeffVals=coeffvalues(fitVals);
predictedGFP=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing))); 
CorrectedGFP=RawData(:,1)-predictedGFP';
FitStruct.GFP=fitVals;
DetrendData=[CorrectedGFP,CorrectedGCaMP,CorrectedCBV];
LowPassData=filtfilt(sos_Low,g_Low,DetrendData);
for q=1:size(LowPassData,2)
RescaleData(:,q)=rescale(LowPassData(:,q),0,1); %rescale all data between 0 to 1
end

%% Model hemodynamic response attenuation of GCaMP signal
BinEdges=(-1:0.005:1);
figure(99);scatter(RescaleData(:,3),RescaleData(:,2));
% figure(99);scatter(LowPassData(:,3),LowPassData(:,2));
% BinEdges=(-1:0.002:1);
figure(100);
DataHist=histogram2(RescaleData(:,3),RescaleData(:,2),BinEdges,BinEdges);
% DataHist=histogram2(LowPassData(:,3),LowPassData(:,2),BinEdges,BinEdges);
BinCounts=DataHist.BinCounts;
ColCounts=sum(BinCounts,2);
NormMat=repmat(ColCounts,1,size(BinCounts,2));
NormHist=(BinCounts./NormMat)*100;
KeepBins=ColCounts>=2*ChunkData.Params.DataFs;


%% Find the bin with the highest number of counts for each real HbT bin
for colNum=1:size(NormHist,1)
    maxVal=max(NormHist(colNum,:));
    if maxVal~=0
        MaxBinLocs=find(NormHist(colNum,:)==maxVal);
        PeakLocs(colNum)=round(median(MaxBinLocs),0);
    else
        PeakLocs(colNum)=NaN;
    end
end
PeakLocs=PeakLocs.*KeepBins';
PeakLocs(PeakLocs==0)=NaN;
PredictedPeaks(1:length(BinEdges))=0;
PredictedPeaks(~isnan(PeakLocs))=BinEdges(PeakLocs(~isnan(PeakLocs)));

figure(103); hold on;
imagesc(BinEdges,BinEdges,NormHist');
caxis([0 15]);
colorbar('eastoutside')
axis xy;
ylabel('Observed GFP');
xlabel('Observed \DeltaCBV');
title('Normalized histogram of CBV vs GFP');
plot(BinEdges,PredictedPeaks,'r','LineWidth',2);
xlim([0 1]);
ylim([0 1]);
% xlim([-0.1 0.3]);
% ylim([-0.1 0.1]);
legend('Peak count at observed \DeltaCBV');
EmptyBins=find(isnan(PeakLocs)==1);
StartPoint=EmptyBins(find(diff(EmptyBins)>5,1,'first'))+1;
if isempty(StartPoint)
    StartPoint=EmptyBins(length(EmptyBins))+1;
end
EndPoint=EmptyBins(find(diff(EmptyBins)>5,1,'first')+1)-1;

if isempty(EndPoint)
    EndPoint=length(ColCounts);
end
[theCurve,goodness,stats]=fit(BinEdges(StartPoint:EndPoint)',PredictedPeaks(StartPoint:EndPoint)','poly1','Robust','LAR');
% [theCurve,goodness,stats]=fit(CenteredData(:,3),CenteredData(:,2),'poly1','Robust','LAR');
coeffVals=coeffvalues(theCurve);
theEqn=formula(theCurve);
theFit=coeffVals(1).*(-0.5:0.001:1)+coeffVals(2);
plot((-0.5:0.001:1),theFit,'c','LineWidth',2);
savefig([ChunkData.Params.animalname '_' ChunkData.Params.date '_NormalizedHistogramFit']);

FinalGCaMP=LowPassData(:,2)-(coeffVals(1)*LowPassData(:,3));
figure(50);plot(LowPassData(:,2));hold on; plot(FinalGCaMP);plot(LowPassData(:,3));
end
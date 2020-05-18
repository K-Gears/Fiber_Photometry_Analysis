function [FigureData]=CreatePlotStruct(GroupedData)
%This function assembles all of the necessary data for plotting figures in
%Echagarugga et al 2020.
%% Useful Constants for data formatting
[z,p,k]=butter(3,0.5/(0.5*GroupedData.CAG_eGFP.Params.DataFs),'low');
[sos_ball,g_ball]=zp2sos(z,p,k);

[z,p,k]=butter(3,1/(0.5*GroupedData.CAG_eGFP.Params.DataFs),'low');
[sos_fbr,g_fbr]=zp2sos(z,p,k);
%% Color Maps
FigureData.ColorMaps.cmap=brewermap(12,'Set2');
FigureData.ColorMaps.AnMap=brewermap(12,'Dark2');
%% Five Second Locomotion Figure
dataTypes=fieldnames(GroupedData);
for dataNum=2:size(dataTypes,1)
PeakVals=squeeze(GroupedData.(dataTypes{dataNum}).AveragedData.five_second_events.PeakResp);  

FigureData.FiveSecFig.ResponseVol.(dataTypes{dataNum}).SingleAnimalAvg.GCaMP=PeakVals(2,:);
FigureData.FiveSecFig.ResponseVol.(dataTypes{dataNum}).SingleAnimalAvg.CBV=PeakVals(3,:);
FigureData.FiveSecFig.ResponseVol.(dataTypes{dataNum}).PopulationAvg.GCaMP=GroupedData.(dataTypes{dataNum}).Averages.AveragedData.five_second_events.PeakResp(2);
FigureData.FiveSecFig.ResponseVol.(dataTypes{dataNum}).PopulationAvg.CBV=GroupedData.(dataTypes{dataNum}).Averages.AveragedData.five_second_events.PeakResp(3);

FigureData.FiveSecFig.LocomotionResponse.(dataTypes{dataNum}).GCaMP=GroupedData.(dataTypes{dataNum}).Averages.AveragedData.five_second_events.AvgResp(:,2);
FigureData.FiveSecFig.LocomotionResponse.(dataTypes{dataNum}).CBV=GroupedData.(dataTypes{dataNum}).Averages.AveragedData.five_second_events.AvgResp(:,3);
FigureData.FiveSecFig.LocomotionResponse.(dataTypes{dataNum}).PlotTime=((1:length(GroupedData.(dataTypes{dataNum}).Averages.AveragedData.five_second_events.AvgResp))-GroupedData.(dataTypes{dataNum}).Params.StartPad)/GroupedData.(dataTypes{dataNum}).Params.DataFs;
end

%% Single Animal Hsyn-GCaMP6s
load('CE_FBR006_Hsyn_GCaMP_1500_um_AnimalAvgData.mat');
load('CE_FBR006_100819_FiberPhotometry.mat');

[~,velocity]=velocity_binarize_fiberphotometry(ChunkData.WheelData.AnalogSignal,ChunkData.Params.DataFs,ChunkData.Params.DataFs,1e-4);
WheelData=filtfilt(sos_ball,g_ball,velocity);%filtfilt(sos_ball,g_ball,abs(diff(velocity)));
% WheelData(WheelData<0)=0;

FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.Trial.PlotTime=(1:length(ChunkData.RawFiberData.ZScoreFiberData))/ChunkData.Params.DataFs;
FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.Trial.BallAcc=WheelData;
FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.Trial.GCaMP=ChunkData.RawFiberData.ZScoreFiberData(:,2);
FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.Trial.CBV=ChunkData.RawFiberData.ZScoreFiberData(:,3);

FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.PlotTime=((1:length(AnimalAvg.AveragedData.five_second_events.AvgResp))-AnimalAvg.Params.StartPad)/AnimalAvg.Params.DataFs;
FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.GCaMP=AnimalAvg.AveragedData.five_second_events.AvgResp(:,2);
FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.GCaMPStd=AnimalAvg.AveragedData.five_second_events.StdResp(:,2);
FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.CBV=AnimalAvg.AveragedData.five_second_events.AvgResp(:,3);
FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.CBVStd=AnimalAvg.AveragedData.five_second_events.StdResp(:,3);

%% Single Animal CaMKII-GCaMP6s
load('CE_FBR022_CaMKII_GCaMP_1500_um_AnimalAvgData.mat');
load('CE_FBR022_030220_FiberPhotometry.mat');

[~,velocity]=velocity_binarize_fiberphotometry(ChunkData.WheelData.AnalogSignal,ChunkData.Params.DataFs,ChunkData.Params.DataFs,1e-4);
 WheelData=filtfilt(sos_ball,g_ball,velocity);%abs(diff(filtfilt(sos_ball,g_ball,velocity)));
% WheelData(WheelData<0)=0;

FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.Trial.PlotTime=(1:length(ChunkData.RawFiberData.ZScoreFiberData))/ChunkData.Params.DataFs;
FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.Trial.BallAcc=WheelData;
FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.Trial.GCaMP=ChunkData.RawFiberData.ZScoreFiberData(:,2);
FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.Trial.CBV=ChunkData.RawFiberData.ZScoreFiberData(:,3);

FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.PlotTime=((1:length(AnimalAvg.AveragedData.five_second_events.AvgResp))-AnimalAvg.Params.StartPad)/AnimalAvg.Params.DataFs;
FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.GCaMP=AnimalAvg.AveragedData.five_second_events.AvgResp(:,2);
FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.GCaMPStd=AnimalAvg.AveragedData.five_second_events.StdResp(:,2);
FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.CBV=AnimalAvg.AveragedData.five_second_events.AvgResp(:,3);
FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.CBVStd=AnimalAvg.AveragedData.five_second_events.StdResp(:,3);

%% Single Animal nNOS-GCaMP6s
load('CE_FBR018_NOS_cre_1500_um_AnimalAvgData.mat');
load('CE_FBR018_013020_FiberPhotometry.mat');

[~,velocity]=velocity_binarize_fiberphotometry(ChunkData.WheelData.AnalogSignal,ChunkData.Params.DataFs,ChunkData.Params.DataFs,1e-4);
WheelData=filtfilt(sos_ball,g_ball,velocity);%abs(diff(filtfilt(sos_ball,g_ball,velocity)));
% WheelData(WheelData<0)=0;

FigureData.FiveSecFig.SingleAnimalPlots.nNOS.Trial.PlotTime=(1:length(ChunkData.RawFiberData.ZScoreFiberData))/ChunkData.Params.DataFs;
FigureData.FiveSecFig.SingleAnimalPlots.nNOS.Trial.BallAcc=WheelData;
FigureData.FiveSecFig.SingleAnimalPlots.nNOS.Trial.GCaMP=ChunkData.RawFiberData.ZScoreFiberData(:,2);
FigureData.FiveSecFig.SingleAnimalPlots.nNOS.Trial.CBV=ChunkData.RawFiberData.ZScoreFiberData(:,3);

FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.PlotTime=((1:length(AnimalAvg.AveragedData.five_second_events.AvgResp))-AnimalAvg.Params.StartPad)/AnimalAvg.Params.DataFs;
FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.GCaMP=AnimalAvg.AveragedData.five_second_events.AvgResp(:,2);
FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.GCaMPStd=AnimalAvg.AveragedData.five_second_events.StdResp(:,2);
FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.CBV=AnimalAvg.AveragedData.five_second_events.AvgResp(:,3);
FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.CBVStd=AnimalAvg.AveragedData.five_second_events.StdResp(:,3);

%% Extended Locomotion Bar Plots
for dataNum=2:size(dataTypes,1)
    FigureData.ExtendedEvents.BarPlots.GCaMPAvg((dataNum-1),:)=[GroupedData.(dataTypes{dataNum}).Averages.AveragedData.ten_second_events.PeakResp(2),GroupedData.(dataTypes{dataNum}).Averages.AveragedData.fifteen_second_events.PeakResp(2),GroupedData.(dataTypes{dataNum}).Averages.AveragedData.thirty_second_events.PeakResp(2)];
    FigureData.ExtendedEvents.BarPlots.CBVAvg((dataNum-1),:)=[GroupedData.(dataTypes{dataNum}).Averages.AveragedData.ten_second_events.PeakResp(3),GroupedData.(dataTypes{dataNum}).Averages.AveragedData.fifteen_second_events.PeakResp(3),GroupedData.(dataTypes{dataNum}).Averages.AveragedData.thirty_second_events.PeakResp(3)];
end 
for dataNum=2:size(dataTypes,1)
    PeakVals10s=squeeze(GroupedData.(dataTypes{dataNum}).AveragedData.ten_second_events.PeakResp);
    PeakVals15s=squeeze(GroupedData.(dataTypes{dataNum}).AveragedData.fifteen_second_events.PeakResp);
    PeakVals30s=squeeze(GroupedData.(dataTypes{dataNum}).AveragedData.thirty_second_events.PeakResp);
    
     FigureData.ExtendedEvents.BarPlots.(dataTypes{dataNum}).GCaMPVals=[PeakVals10s(2,:);PeakVals15s(2,:);PeakVals30s(2,:)];
     FigureData.ExtendedEvents.BarPlots.(dataTypes{dataNum}).CBVVals=[PeakVals10s(3,:);PeakVals15s(3,:);PeakVals30s(3,:)];
end

%% Ten Second Event Averages
for dataNum=2:size(dataTypes,1)
    FigureData.ExtendedEvents.TenSecEventAvg.(dataTypes{dataNum}).PlotTime=((1:length(GroupedData.(dataTypes{dataNum}).Averages.AveragedData.ten_second_events.AvgResp))-GroupedData.(dataTypes{dataNum}).Params.StartPad)/GroupedData.(dataTypes{dataNum}).Params.DataFs;
    FigureData.ExtendedEvents.TenSecEventAvg.(dataTypes{dataNum}).GCaMP=GroupedData.(dataTypes{dataNum}).Averages.AveragedData.ten_second_events.AvgResp(:,2);
    FigureData.ExtendedEvents.TenSecEventAvg.(dataTypes{dataNum}).CBV=GroupedData.(dataTypes{dataNum}).Averages.AveragedData.ten_second_events.AvgResp(:,3);
end

%% Fifteen Second Event Averages
for dataNum=2:size(dataTypes,1)
    FigureData.ExtendedEvents.FifteenSecEventAvg.(dataTypes{dataNum}).PlotTime=((1:length(GroupedData.(dataTypes{dataNum}).Averages.AveragedData.fifteen_second_events.AvgResp))-GroupedData.(dataTypes{dataNum}).Params.StartPad)/GroupedData.(dataTypes{dataNum}).Params.DataFs;
    FigureData.ExtendedEvents.FifteenSecEventAvg.(dataTypes{dataNum}).GCaMP=GroupedData.(dataTypes{dataNum}).Averages.AveragedData.fifteen_second_events.AvgResp(:,2);
    FigureData.ExtendedEvents.FifteenSecEventAvg.(dataTypes{dataNum}).CBV=GroupedData.(dataTypes{dataNum}).Averages.AveragedData.fifteen_second_events.AvgResp(:,3);
end

%% Thirty Second Event Averages
for dataNum=2:size(dataTypes,1)
    FigureData.ExtendedEvents.ThirtySecEventAvg.(dataTypes{dataNum}).PlotTime=((1:length(GroupedData.(dataTypes{dataNum}).Averages.AveragedData.thirty_second_events.AvgResp))-GroupedData.(dataTypes{dataNum}).Params.StartPad)/GroupedData.(dataTypes{dataNum}).Params.DataFs;
    FigureData.ExtendedEvents.ThirtySecEventAvg.(dataTypes{dataNum}).GCaMP=GroupedData.(dataTypes{dataNum}).Averages.AveragedData.thirty_second_events.AvgResp(:,2);
    FigureData.ExtendedEvents.ThirtySecEventAvg.(dataTypes{dataNum}).CBV=GroupedData.(dataTypes{dataNum}).Averages.AveragedData.thirty_second_events.AvgResp(:,3);
end

%% Cross Correlation
for dataNum=2:size(dataTypes,1)
    FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).Lags=GroupedData.(dataTypes{dataNum}).Averages.Xcorr.Lags;
    FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).CorrCoeff=GroupedData.(dataTypes{dataNum}).Averages.Xcorr.CBV_GCaMP;
    
for anNum=1:size(GroupedData.(dataTypes{dataNum}).Xcorr.CBV_GCaMP,3)
        PeakLoc=find(GroupedData.(dataTypes{dataNum}).Xcorr.CBV_GCaMP(:,1,anNum)==max(GroupedData.(dataTypes{dataNum}).Xcorr.CBV_GCaMP(:,1,anNum),[],'all'));
        tempInds=GroupedData.(dataTypes{dataNum}).Xcorr.Lags(:,:,anNum)>=-5 & GroupedData.(dataTypes{dataNum}).Xcorr.Lags(:,:,anNum)<=5;
        PeakAlt=find(diff(GroupedData.(dataTypes{dataNum}).Xcorr.CBV_GCaMP(tempInds,1,anNum),1,1)<0,1,'first')+(find(tempInds==1,1,'first')-1);
        if PeakAlt<PeakLoc
            PeakLoc=PeakAlt;
        end

    PeakTime(anNum)=GroupedData.(dataTypes{dataNum}).Xcorr.Lags(1,PeakLoc,anNum);
    PeakVal(anNum)=GroupedData.(dataTypes{dataNum}).Xcorr.CBV_GCaMP(PeakLoc,1,anNum);
end

FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).PeakTime=PeakTime;
FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).PeakVal=PeakVal;
FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).AvgTime=mean(PeakTime);
FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).AvgVal=mean(PeakVal);
clear PeakTime PeakVal
end

%% Coherence
for dataNum=2:size(dataTypes,1)
    FigureData.Xcorr.Coherence.(dataTypes{dataNum}).Freqs=GroupedData.(dataTypes{dataNum}).Averages.FrequencyDomain.Frequency_BloodVolume_GCaMP6s;
    FigureData.Xcorr.Coherence.(dataTypes{dataNum}).Coherence=GroupedData.(dataTypes{dataNum}).Averages.FrequencyDomain.Coherence_BloodVolume_GCaMP6s;
end

%% Signal to Noise
    for dataNum=2:size(dataTypes,1)
        FigureData.Xcorr.SNR.GCaMPAvg((dataNum-1),:)=GroupedData.(dataTypes{dataNum}).Averages.RawFiberData.RawSignalAvg(2);
        FigureData.Xcorr.SNR.GCaMPSTD((dataNum-1),:)=GroupedData.(dataTypes{dataNum}).Averages.RawFiberData.LowPassSignalSTD(2);
        FigureData.Xcorr.SNR.CBVAvg((dataNum-1),:)=GroupedData.(dataTypes{dataNum}).Averages.RawFiberData.RawSignalAvg(3);
        FigureData.Xcorr.SNR.CBVSTD((dataNum-1),:)=GroupedData.(dataTypes{dataNum}).Averages.RawFiberData.LowPassSignalSTD(3);
    end
    
    for dataNum=2:size(dataTypes,1)
        RawSignalAvg=squeeze(GroupedData.(dataTypes{dataNum}).RawFiberData.RawSignalAvg);
        LowPassSTD=squeeze(GroupedData.(dataTypes{dataNum}).RawFiberData.LowPassSignalSTD);
        FigureData.Xcorr.SNR.(dataTypes{dataNum}).RawGCaMPAvg= RawSignalAvg(2,:);
        FigureData.Xcorr.SNR.(dataTypes{dataNum}).RawGCaMPSTD=LowPassSTD(2,:);
        FigureData.Xcorr.SNR.(dataTypes{dataNum}).RawCBVAvg=RawSignalAvg(3,:);
        FigureData.Xcorr.SNR.(dataTypes{dataNum}).RawCBVSTD=LowPassSTD(3,:);
        FigureData.Xcorr.SNR.(dataTypes{dataNum}).GenotypeRawPDF=normpdf((0:0.001:7),GroupedData.(dataTypes{dataNum}).Averages.RawFiberData.RawSignalAvg(2),GroupedData.(dataTypes{dataNum}).Averages.RawFiberData.LowPassSignalSTD(2));
        FigureData.Xcorr.SNR.(dataTypes{dataNum}).GenotypeRawNormPDF=FigureData.Xcorr.SNR.(dataTypes{dataNum}).GenotypeRawPDF/max(FigureData.Xcorr.SNR.(dataTypes{dataNum}).GenotypeRawPDF);
    end

 %% Signal Correction
load('CE_FBR001_CAG_eGFP_1500_um_AnimalAvgData.mat');
load('CE_FBR001_082119_FiberPhotometry.mat');

%% Correct CBV
coeffValsExp=coeffvalues(FitStruct.CBV);
Spacing=1:1:length(ChunkData.RawFiberData.RawFiberData(:,2));
predictedCBVExp=(coeffValsExp(1)*exp((coeffValsExp(2).*Spacing)))+(coeffValsExp(3)*exp((coeffValsExp(4).*Spacing)));
CorrectedCBV=ChunkData.RawFiberData.RawFiberData(:,3)-predictedCBVExp';

FigureData.SignalCorrection.CorrectMethod.PlotTime=(1:length(ChunkData.RawFiberData.RawFiberData))/ChunkData.Params.DataFs;
FigureData.SignalCorrection.CorrectMethod.RawCBV=ChunkData.RawFiberData.RawFiberData(:,3);
FigureData.SignalCorrection.CorrectMethod.PredictedCBVDecay=(coeffValsExp(1)*exp((coeffValsExp(2).*Spacing)))+(coeffValsExp(3)*exp((coeffValsExp(4).*Spacing)));
FigureData.SignalCorrection.CorrectMethod.CorrectedCBV=filtfilt(sos_fbr,g_fbr,CorrectedCBV);

%% Correct GCaMP
coeffValsExp=coeffvalues(FitStruct.GCaMP);
Spacing=1:1:length(ChunkData.RawFiberData.RawFiberData(:,2));
predictedGCaMPExp=(coeffValsExp(1)*exp((coeffValsExp(2).*Spacing)))+(coeffValsExp(3)*exp((coeffValsExp(4).*Spacing)));
CorrectedGCaMP=ChunkData.RawFiberData.RawFiberData(:,2)-predictedGCaMPExp';

FigureData.SignalCorrection.CorrectMethod.RawGCaMP=ChunkData.RawFiberData.RawFiberData(:,2);
FigureData.SignalCorrection.CorrectMethod.PredictedGCaMPDecay=(coeffValsExp(1)*exp((coeffValsExp(2).*Spacing)))+(coeffValsExp(3)*exp((coeffValsExp(4).*Spacing)));
FigureData.SignalCorrection.CorrectMethod.CorrectedGCaMP=filtfilt(sos_fbr,g_fbr,CorrectedGCaMP);

%% Single Animal Correction eGFP
load('CE_FBR001_CAG_eGFP_1500_um_AnimalAvgData.mat');
load('CE_FBR001_082119_FiberPhotometry.mat');
[~,velocity]=velocity_binarize_fiberphotometry(ChunkData.WheelData.AnalogSignal,ChunkData.Params.DataFs,ChunkData.Params.DataFs,1e-4);
WheelData=filtfilt(sos_ball,g_ball,(velocity+8.434e-4));%abs(diff(filtfilt(sos_ball,g_ball,velocity)));
WheelData(WheelData<0)=0;

FigureData.SignalCorrection.eGFPTrial.PlotTime=(1:length(ChunkData.RawFiberData.ZScoreFiberData))/ChunkData.Params.DataFs;
FigureData.SignalCorrection.eGFPTrial.BallAcc=WheelData;
FigureData.SignalCorrection.eGFPTrial.UncorrectedGCaMP=ChunkData.RawFiberData.UncorrectedZScoreGCaMPData((60*ChunkData.Params.DataFs):end);
FigureData.SignalCorrection.eGFPTrial.CorrectedGCaMP=ChunkData.RawFiberData.ZScoreFiberData(:,2);
FigureData.SignalCorrection.eGFPTrial.CorrectedCBV=ChunkData.RawFiberData.ZScoreFiberData(:,3);

%% Correction Histogram
load('CE_FBR001_082119_FiberPhotometry.mat');
[BinEdges,NormHist,theFit]=Correction_Histogram(ChunkData);
FigureData.SignalCorrection.CorrectionHist.BinEdges=BinEdges;
FigureData.SignalCorrection.CorrectionHist.BinCounts=NormHist;
FigureData.SignalCorrection.CorrectionHist.LinearFit=theFit;
end

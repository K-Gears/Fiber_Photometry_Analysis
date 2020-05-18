function [ChunkData]= ExtractFPdata_003(filename,OpticalChannelNames,AnalogChannelNames,CorrectionConst)
%% READ ME
%function used to extract imaging data from .csv files generated during
%fiber photometry experiments
%filename: string of filename to be read
%OpticalChannelNames:Cell array containing names to be assigned to optical
%channels. {'autofluoresence','EGFP','TRITC'}
%AnalogChannelNames:Cell array containing names to be assigned to analog
%signals {'wheel','left whisker solenoid','right whisker solenoid}

close all
%% Constants
ChunkData.Params.Channels=[OpticalChannelNames,AnalogChannelNames];
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
ChunkData.Params.Fit_Freq=0.05;
ChunkData.Params.low_Freq=1;
ChunkData.Params.Final_Freq=[0.01 1];
ChunkData.Params.Correct_Freq=0.01;
params.Fs=ChunkData.Params.DataFs; %multitaper estimation parameters
params.tapers=[3 5];%multitaper estimation parameters

[z,p,k]=butter(3,ChunkData.Params.Fit_Freq/(0.5*ChunkData.Params.DataFs),'low'); %design lowpass filter for hemodynamic correction
[sos_Fit,g_Fit]=zp2sos(z,p,k);

[z,p,k]=butter(3,ChunkData.Params.Correct_Freq/(0.5*ChunkData.Params.DataFs),'low'); %design lowpass filter for hemodynamic correction
[sos_Correct,g_Correct]=zp2sos(z,p,k);

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
[imp_bin]=velocity_binarize(WheelData,ChunkData.Params.DataFs,ChunkData.Params.DataFs,1e-3);
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
[fitVals]=fit(Spacing',RawData(:,3),'exp2','Exclude',ExcludeVals);
coeffVals=coeffvalues(fitVals);
predictedCBV=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing))); 
% figTime=(1:length(RawData(:,3)))/ChunkData.Params.DataFs;
% figure(1);plot(figTime,RawData(:,3)); hold on; plot(figTime,predictedCBV); plot(figTime,FiltData(:,3));title('Exponential fit of TRITC metabolism'); xlabel('Time (sec)'); legend({'Raw TRITC brightness','Exponential Fit','Low pass filtered data fit'}); xlim([0 figTime(end)]);
% xticks(1:900:figTime(length(figTime)));
CorrectedCBV=RawData(:,3)-predictedCBV';
FitStruct.CBV=fitVals;
% Correct Ca2+ dependent GCaMP
[fitVals]=fit(Spacing',RawData(:,2),'exp2','Exclude',ExcludeVals);
coeffVals=coeffvalues(fitVals);
predictedGCaMP=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing))); 
% figure(5);plot(figTime,RawData(:,2)); hold on; plot(figTime,predictedGCaMP);plot(figTime,FiltData(:,2)); title('Exponential fit of GCaMP6s photobleaching'); xlabel('Time (sec)'); legend({'Raw GCaMP brightness','Exponential Fit','Low pass filtered data fit'}); xlim([0 figTime(end)]);
% xticks(1:900:figTime(length(figTime)));
CorrectedGCaMP=RawData(:,2)-predictedGCaMP';
FitStruct.GCaMP=fitVals;
% Correct Ca2+ independent GCaMP
[fitVals]=fit(Spacing',RawData(:,1),'exp2','Exclude',ExcludeVals);
coeffVals=coeffvalues(fitVals);
predictedGFP=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing))); 
% figure(6);plot(figTime,RawData(:,1)); hold on; plot(figTime,predictedGFP); plot(figTime,FiltData(:,1)); title('Exponential fit of GCaMP6s photobleaching'); xlabel('Time (sec)'); legend({'Raw GCaMP brightness','Exponential Fit','Low pass filtered data fit'}); xlim([0 figTime(end)]);
% xticks(1:900:figTime(length(figTime)));
CorrectedGFP=RawData(:,1)-predictedGFP';
FitStruct.GFP=fitVals;
DetrendData=[CorrectedGFP,CorrectedGCaMP,CorrectedCBV];
LowPassData=filtfilt(sos_Low,g_Low,DetrendData);

%% Model hemodynamic response attenuation of GCaMP signal
% for num=1:size(DetrendData,2)
% FitPassData(:,num)=DetrendData(:,num);
% end
% %filtfilt(sos_Fit,g_Fit,DetrendData); %filter data below 0.05Hz to fit hemodynamics to GFP control signal
% % [fitVals]=fit(FitPassData(:,3),FitPassData(:,1),'poly2');
% % figure(98); plot(fitVals,FitPassData(:,3),FitPassData(:,1));
% % varNames=coeffnames(fitVals);
% % Eqn=formula(fitVals);
% % coeffValsPoly=coeffvalues(fitVals);
% 
% [fitVals]=fit(FitPassData(:,3),FitPassData(:,1),'exp2');
% figure(99); plot(fitVals,FitPassData(:,3),FitPassData(:,1));
% varNames=coeffnames(fitVals);
% Eqn=formula(fitVals);
% coeffValsExp=coeffvalues(fitVals);
% FitStruct.Hemo=fitVals;
% predictedGFPExp=(coeffValsExp(1)*exp((coeffValsExp(2).*LowPassData(:,3))))+(coeffValsExp(3)*exp((coeffValsExp(4).*LowPassData(:,3)))); %exponential fit
% ScalePredictedGFP=rescale(predictedGFPExp,-1,1);
%predictedGFPPoly=coeffValsPoly(1)*DetrendData(:,3).^2+coeffValsPoly(2)*DetrendData(:,3)+coeffValsPoly(3); %polynomial fit

%ZscorePredictedPoly=filtfilt(sos_Fit,g_Fit,(predictedGFPPoly-mean(predictedGFPPoly))/std(predictedGFPPoly));
% ZscorePredictedExp=filtfilt(sos_Fit,g_Fit,(predictedGFPExp-mean(predictedGFPExp))/std(predictedGFPExp));
% ZscorePredictedExp=(predictedGFPExp-mean(predictedGFPExp))/std(predictedGFPExp);
%RsqrPredicted=1-(sum((LowPassData(:,1)-predictedGFP).^2)/sum((LowPassData(:,1)-mean(LowPassData(:,1))).^2));

% figure(2);plot(figTime,FitPassData(:,1));hold on; plot(figTime,predictedGFPExp); title('Exponential fit of GFP brightness vs CBV'); xlabel('Time (sec)'); legend({'Raw GFP brightness','Exponential Fit'}); xlim([0 figTime(end)]);%plot(figTime,predictedGFPPoly);

%% Correct for Hemodynamic attenuation by subtracting 560nm signal from 465nm
OffsetStart=15*ChunkData.Params.DataFs; %Removes filtering artifact at start of "Smooth465"
% FiltData=filtfilt(sos_Low,g_Low,DetrendData); %Low pass filter all data below 1Hz
for q=1:size(LowPassData,2)
RescaleData(:,q)=rescale(LowPassData(:,q),0,1); %rescale all data between 0 to 1
% RescaleData(:,q)=Tempscale-mean(Tempscale);
end
% Corrected465=RescaleData(:,2)-(CorrectionConst*RescaleData(:,3)); % Add rescaled TRITC signal to GCaMP signal
Corrected465=RescaleData(:,2)-(CorrectionConst*RescaleData(:,3)); % Add rescaled TRITC signal to GCaMP signal
Smooth465=filtfilt(sos_final,g_final,Corrected465); %Bandpass filter data between [0.01 and 1] Hz
Z465=(Smooth465-mean(Smooth465(OffsetStart:end)))/std(Smooth465(OffsetStart:end)); % Z score GCaMP data
% SmoothData=filtfilt(sos_final,g_final,RescaleData);%Bandpass filter data between [0.01 and 1] Hz
SmoothData=filtfilt(sos_final,g_final,RescaleData);%Bandpass filter data between [0.01 and 1] Hz
Z560=(SmoothData(:,3)-mean(SmoothData(:,3)))/std(SmoothData(:,3)); % Z score TRITC data


%% Z-score optical data
% AvgData=mean(RescaleData,1);
% StdData=std(RescaleData,0,1);
AvgData=mean(RescaleData,1);
StdData=std(RescaleData,0,1);
AvgMatrix=repmat(AvgData,length(RescaleData),1);
StdMatrix=repmat(StdData,length(RescaleData),1);
% ZscoredFiberData=(RescaleData-AvgMatrix)./StdMatrix;
ZscoredFiberData=(RescaleData-AvgMatrix)./StdMatrix;
% figure(7);plot(figTime,ZscoredFiberData(:,3));hold on;  plot(figTime,ZscorePredictedExp);title('Exponential fit of GFP brightness vs CBV'); xlabel('Time (sec)'); legend({'Z-score TRITC brightness','Z-score GFP Exponential Fit'}); xlim([0 figTime(end)]);%plot(figTime,ZscorePredictedPoly);

%% Correct GCaMP channel for hemodynamic attenuation

% figure(3);plot(figTime,ZscoredFiberData(:,2)); hold on; plot(figTime,ZscorePredictedExp); legend({'Raw GCaMP Z Score','Exponential fit hemodynamic attenuation'});title('Predicted hemodynamic evoked GCaMP6s attenuation and un-corrected GCaMP6s Z-score'); ylabel('Z-score');xlabel('Time (sec)');xlim([0 figTime(end)]);
% RawGCaMPZ=ZscoredFiberData(:,2);
% ZscoredFiberData(:,2)=ZscoredFiberData(:,2)-ZscorePredictedExp;
% correctedGCaMPz=filtfilt(sos_final,g_final,ZscoredFiberData(:,2));%ZscoredFiberData(:,2);
% figure(4);hold on;plot(figTime,ZscoredFiberData(:,2)); plot(figTime,correctedGCaMPz);plot(figTime,Z465);  legend({'Raw GCaMP Z-score','Hemodynamic corrected GCaMP6s Z-score','Bandpass filtered hemodynamic corrected','560nm subtraction Correction'});title('Correction of hemodynamic evoked GCaMP6s attenuation'); ylabel('Z-score');xlabel('Time (sec)');xlim([0 figTime(end)]);
% difference465=Z465-correctedGCaMPz;
UncorrectedZscoredGCaMPData=filtfilt(sos_final,g_final,ZscoredFiberData(:,2));    
ZscoredFiberData(:,2)=Z465;%filtfilt(sos_final,g_final,ZscoredFiberData(:,2));
ZscoredFiberData(:,3)=filtfilt(sos_Low,g_Low,ZscoredFiberData(:,3));
% figure(101);plot(figTime,ZscoredFiberData(:,(2:3)));hold on; plot(figTime,correctedGCaMPz);xlabel('Time (sec)'); ylabel('z-score'); legend({'560nm subtracted 465nm GCaMP','Blood volume','405nm Corrected 465nm GCaMP'}); title('Trial GCaMP and CBV');
% difference560=Z560-ZscoredFiberData(:,3);
% figure;plot(figTime,difference465);
% figure;plot(difference560);



%% Trim to eliminate filtering ring
StartInd=60*ChunkData.Params.DataFs;
WheelData=WheelData(StartInd:end);
ZscoredFiberData=ZscoredFiberData((StartInd:end),:);
LowPassData=LowPassData((StartInd:end),:);
%% Cross Spectral Frequency Coherence
params.fpass=[0.01 0.5];
params.tapers=[19 37];

[C,phi,S12,S1,S2,f]=coherencyc(Z560,Z465,params);
ChunkData.FrequencyDomain.(['Coherence_' OpticalChannelNames{3} '_' OpticalChannelNames{2}])=C;
ChunkData.FrequencyDomain.(['Phase_' OpticalChannelNames{3} '_' OpticalChannelNames{2}])=phi;
ChunkData.FrequencyDomain.(['Frequency_' OpticalChannelNames{3} '_' OpticalChannelNames{2}])=f;

% [C,phi,S12,S1,S2,f]=coherencyc(ZscoredFiberData(:,3),ZscoredFiberData(:,1),params);
% ChunkData.FrequencyDomain.(['Coherence_' OpticalChannelNames{3} '_' OpticalChannelNames{1}])=C;
% ChunkData.FrequencyDomain.(['Phase_' OpticalChannelNames{3} '_' OpticalChannelNames{1}])=phi;
% ChunkData.FrequencyDomain.(['Frequency_' OpticalChannelNames{3} '_' OpticalChannelNames{1}])=f;

% [C,phi,S12,S1,S2,f]=coherencyc(ZscoredFiberData(:,2),ZscoredFiberData(:,1),params);
% ChunkData.FrequencyDomain.(['Coherence_' OpticalChannelNames{2} '_' OpticalChannelNames{1}])=C;
% ChunkData.FrequencyDomain.(['Phase_' OpticalChannelNames{2} '_' OpticalChannelNames{1}])=phi;
% ChunkData.FrequencyDomain.(['Frequency_' OpticalChannelNames{2} '_' OpticalChannelNames{1}])=f;

%% Binarize locomotion data
[imp_bin]=velocity_binarize(WheelData,ChunkData.Params.DataFs,ChunkData.Params.DataFs,1e-4);

%% Cross correlation analysis
ExpectedLength=ChunkData.Params.DataFs*ChunkData.Params.DataSeconds;
maxLag=60*ChunkData.Params.DataFs; %Calculate cross correlation with +/- 5s lags

[r_GCaMP_Loco]=xcorr(detrend(ZscoredFiberData((1:ExpectedLength),2),'linear'),detrend(double(imp_bin(1:ExpectedLength)),'linear'),maxLag,'coeff');

[r_CBV_Loco]=xcorr(detrend(ZscoredFiberData((1:ExpectedLength),3),'linear'),detrend(double(imp_bin(1:ExpectedLength)),'linear'),maxLag,'coeff');

% [r_CBV_GCaMP]=xcorr(detrend(ZscoredFiberData((1:ExpectedLength),3),'linear'),detrend(ZscoredFiberData((1:ExpectedLength),2),'linear'),maxLag,'coeff');

[r_CBV_GCaMP,lags]=xcorr(detrend(Z560,'linear'),detrend(Z465,'linear'),maxLag,'coeff');
% [r_CBV_GFP,lags]=xcorr(detrend(ZscoredFiberData((1:ExpectedLength),3),'linear'),detrend(ZscoredFiberData((1:ExpectedLength),1),'linear'),maxLag,'coeff');

%% Run triggered averaging
RunDuration=[5 10 15 30 45];
RunLengths={'five_second_events','ten_second_events','fifteen_second_events','thirty_second_events','fortyfive_second_events'};
LeadTime=5;    
[~,~,new_T_run]=motion_cont_2(imp_bin,ChunkData.Params.DataFs,RunDuration(1),LeadTime);
EventLengths=(new_T_run(2,:)-new_T_run(1,:))/ChunkData.Params.DataFs;
for durNum=1:length(RunDuration)
    LocomotionEvents=[];
    RunLabel=RunLengths{durNum};
    
    if durNum==length(RunDuration)
        EventFind=EventLengths>=RunDuration(durNum);
        EventInds=new_T_run(:,EventFind);
    else
        EventFind=EventLengths>=RunDuration(durNum) & EventLengths<RunDuration(durNum+1);
        EventInds=new_T_run(:,EventFind);
    end
    
    for k=1:size(EventInds,2)
        StartInd=EventInds(1,k)-ChunkData.Params.StartPad;
        if durNum==length(RunDuration)
            EndInd=EventInds(1,k)+(ChunkData.Params.DataFs*60);
        else
            EndInd=EventInds(1,k)+(ChunkData.Params.DataFs*RunDuration(durNum+1));%ChunkData.Params.FollowPad;
        end
        if StartInd>0
            if EndInd<length(ZscoredFiberData)+1
                Baseline=mean(ZscoredFiberData((StartInd:(StartInd+(2*ChunkData.Params.DataFs))),:),1);
                BaselineRaw=mean(LowPassData((StartInd:(StartInd+(2*ChunkData.Params.DataFs))),:),1);
                for j=1:size(Baseline,2)
                    LocomotionEvents(:,j,k)=ZscoredFiberData((StartInd:EndInd),j)-Baseline(j);
                    LocomotionEventsRaw(:,j,k)=LowPassData((StartInd:EndInd),j)-BaselineRaw(j);
                end
            end
        end
    end
    
    if ~isempty(LocomotionEvents)
        AverageLocomotionResponse=mean(LocomotionEvents,3);
        StanDevLocomotion=std(LocomotionEvents,0,3);
        MedianLocomotionResponse=median(LocomotionEvents,3);
        BaseResp=mean(AverageLocomotionResponse((1:(3*ChunkData.Params.DataFs)),:),1);
        PeakVal=max(AverageLocomotionResponse,[],1);
        for chanNum=1:size(AverageLocomotionResponse,2)
            PeakTime(chanNum)=(find(AverageLocomotionResponse(:,chanNum)==PeakVal(chanNum))-ChunkData.Params.StartPad)/ChunkData.Params.DataFs;
        end
        PeakInds=((ChunkData.Params.StartPad+RunDuration(durNum)*ChunkData.Params.DataFs)-(3*ChunkData.Params.DataFs)):(ChunkData.Params.StartPad+RunDuration(durNum)*ChunkData.Params.DataFs);
        PeakResp=mean(AverageLocomotionResponse(PeakInds,:),1);
        
        ChunkData.WheelData.(RunLabel).RunInds=EventInds;
        ChunkData.WheelData.(RunLabel).RunLengths=(EventInds(2,:)-EventInds(1,:))/ChunkData.Params.DataFs;
        ChunkData.AveragedData.(RunLabel).AvgResp=AverageLocomotionResponse;
        ChunkData.AveragedData.(RunLabel).StdResp=StanDevLocomotion;
        ChunkData.AveragedData.(RunLabel).MedResp=MedianLocomotionResponse;
        ChunkData.AveragedData.(RunLabel).BaseResp=BaseResp;
        ChunkData.AveragedData.(RunLabel).PeakResp=PeakResp;
        ChunkData.AveragedData.(RunLabel).PeakVal=PeakVal;
        ChunkData.AveragedData.(RunLabel).PeakTime=PeakTime;
        ChunkData.LocomotionEvokedData.(RunLabel).OpticalData=LocomotionEvents;
        ChunkData.LocomotionEvokedData.(RunLabel).OpticalDataRaw=LocomotionEventsRaw;
    else
        ChunkData.WheelData.(RunLabel).RunInds=[];
        ChunkData.WheelData.(RunLabel).RunLengths=[];
        ChunkData.AveragedData.(RunLabel).AvgResp=[];
        ChunkData.AveragedData.(RunLabel).StdResp=[];
        ChunkData.AveragedData.(RunLabel).MedResp=[];
        ChunkData.AveragedData.(RunLabel).BaseResp=[];
        ChunkData.AveragedData.(RunLabel).PeakResp=[];
        ChunkData.AveragedData.(RunLabel).PeakVal=[];
        ChunkData.AveragedData.(RunLabel).PeakTime=[];
        ChunkData.LocomotionEvokedData.(RunLabel).OpticalData=[];
    end
    
%     plotTime=((1:length(ChunkData.AveragedData.(RunLabel).AvgResp))-ChunkData.Params.StartPad)/ChunkData.Params.DataFs;
%     figure;plot(plotTime,ChunkData.AveragedData.(RunLabel).AvgResp(:,(2:3))); xlim([-5 RunDuration(durNum)]); legend({'nNosCre GCaMP6s','TRITC Blood Volume'});
%     title('Average locomotion evoked response'); xlabel('Time (sec)'); ylabel('Z score');
    clear LocomotionEvents LocomotionEventsRaw
end
ChunkData.RawFiberData.ZScoreFiberData=ZscoredFiberData;
ChunkData.RawFiberData.UncorrectedZScoreGCaMPData=UncorrectedZscoredGCaMPData;
ChunkData.RawFiberData.DetrendData=DetrendData;
ChunkData.RawFiberData.RawFiberData=RawData;
ChunkData.RawFiberData.LowPassData=LowPassData;
ChunkData.RawFiberData.RescaledData=RescaleData;
ChunkData.RawFiberData.Corrected465=Corrected465;
ChunkData.WheelData.AnalogSignal=WheelData;
ChunkData.WheelData.BinarizedLocomotion=imp_bin;
ChunkData.Xcorr.Loco_GCaMP=r_GCaMP_Loco;
ChunkData.Xcorr.Loco_CBV=r_CBV_Loco;
ChunkData.Xcorr.CBV_GCaMP=r_CBV_GCaMP;
% ChunkData.Xcorr.CBV_GFP=r_CBV_GFP;
ChunkData.Xcorr.Lags=lags/ChunkData.Params.DataFs;

save([ChunkData.Params.animalname '_' ChunkData.Params.date '_FiberPhotometry.mat'],'ChunkData','FitStruct','-v7.3');


% plotTime=((1:length(ChunkData.AveragedData.MedResp))-ChunkData.Params.StartPad)/ChunkData.Params.DataFs;
% figure;plot(plotTime,ChunkData.AveragedData.MedResp(:,(2:3))); xlim([-5 15]); legend({'nNosCre GCaMP6s','TRITC Blood Volume'});title('Median Locomotion evoked response');xlabel('Time (sec)');ylabel('Z score');

% figure;plot(ChunkData.Xcorr.Lags,ChunkData.Xcorr.CBV_GCaMP);xlim([-10 10]); title('Cross-correlation blood volume x GCaMP6s'); xlabel('Time (sec)'); ylabel('Normalized corr coeff');
% 
% runPlot(1:length(ChunkData.RawFiberData.ZScoreFiberData(:,3)))=0;
% for runNum=1:length(ChunkData.WheelData.five_second_events.RunInds)
%     runPlot(ChunkData.WheelData.five_second_events.RunInds(1,runNum):ChunkData.WheelData.five_second_events.RunInds(2,runNum))=1.1*max(ChunkData.RawFiberData.ZScoreFiberData(:,3));
% end
% runPlot(runPlot==0)=NaN;
% plotTime=(1:length(ChunkData.RawFiberData.ZScoreFiberData(:,3)))/ChunkData.Params.DataFs;
% 
% figure;plot(plotTime,ChunkData.RawFiberData.ZScoreFiberData(:,(2:3)));title('Signal intensity changes');xlabel('Time (sec)');ylabel('Z-score'); 
% xlim([0 plotTime(end)]);
% hold on; scatter(plotTime,runPlot,'k','filled'); legend({'Ca2+ GCaMP6s','TRITC Blood volume','Locomotion events'});
end



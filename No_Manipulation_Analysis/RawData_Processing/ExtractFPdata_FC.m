function [ChunkData]= ExtractFPdata_FC(filename,OpticalChannelNames,AnalogChannelNames,trackingData,correctFlag)
%% READ ME
%function used to extract imaging data from .csv files generated during
%fiber photometry experiments
%filename: string of filename to be read
%OpticalChannelNames:Cell array containing names to be assigned to optical
%channels. {'autofluoresence','EGFP','TRITC'}
%AnalogChannelNames:Cell array containing names to be assigned to analog
%signals {'wheel','left whisker solenoid','right whisker solenoid}
% trackingData: structure containing ROI based movement traces and event
% indicies defined in 'processDoricBehaviorCamera.m'

close all
%% Read .CSV
fileInf=detectImportOptions(filename);
DataChannels=[2,3,5];%these are the three demodulated optical channels [2,3,5] and the wheel data [7]
RawData=csvread(filename,2,0);


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
ChunkData.Params.Acquisition_Fs=1.2e4; %This is evidently highly error prone and noisy
ChunkData.Params.Decimation=10;
ChunkData.Params.DataAcquired=2;%time in hours of data acquired
ChunkData.Params.DataSeconds=ChunkData.Params.DataAcquired*(60^2);
ChunkData.Params.DataFs=round(length(RawData)/ChunkData.Params.DataSeconds,0);%ChunkData.Params.Acquisition_Fs/ChunkData.Params.Decimation;
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

%% Get ball movement and optical channel data
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
[fitVals]=fit(Spacing',RawData(:,3),'exp2','Exclude',ExcludeVals);
coeffVals=coeffvalues(fitVals);
predictedCBV=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing)));
figTime=(1:length(RawData(:,3)))/ChunkData.Params.DataFs;
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

%% Correct for Hemodynamic attenuation by subtracting 560nm signal from 465nm
OffsetStart=60*ChunkData.Params.DataFs; %Removes filtering artifact at start of "Smooth465"
% FiltData=filtfilt(sos_Low,g_Low,DetrendData); %Low pass filter all data below 1Hz
for q=1:size(LowPassData,2)
    RescaleData(:,q)=rescale(LowPassData(:,q),0,1); %rescale all data between 0 to 1
end
if strcmpi(correctFlag,'y')
    [correctionConstant,initCorrCoeffs,adjCorrCoeffs]= MinimizeCorrCoeff(RescaleData(:,2),RescaleData(:,3),ChunkData.Params.DataFs);
    Corrected465=RescaleData(:,2)+(correctionConstant*RescaleData(:,3));
else
    Corrected465=RescaleData(:,2);
end
% Corrected465=RescaleData(:,2)-(CorrectionConst*RescaleData(:,3)); % Add rescaled TRITC signal to GCaMP signal
Smooth465=filtfilt(sos_Low,g_Low,Corrected465);%RescaleData(:,2)); %Bandpass filter data between [0.01 and 1] Hz
Z465=(Smooth465-mean(Smooth465(OffsetStart:end)))/std(Smooth465(OffsetStart:end)); % Z score GCaMP data
% SmoothData=filtfilt(sos_final,g_final,RescaleData);%Bandpass filter data between [0.01 and 1] Hz
SmoothData=filtfilt(sos_Low,g_Low,RescaleData);%Bandpass filter data between [0.01 and 1] Hz
Z560=(SmoothData(:,3)-mean(SmoothData((OffsetStart:end),3)))/std(SmoothData((OffsetStart:end),3)); % Z score TRITC data

%% Z-score optical data
AvgData=mean(RescaleData,1);
StdData=std(RescaleData,0,1);
AvgMatrix=repmat(AvgData,length(RescaleData),1);
StdMatrix=repmat(StdData,length(RescaleData),1);
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
ZscoredFiberData(:,3)=Z560;%filtfilt(sos_Low,g_Low,ZscoredFiberData(:,3));
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

% [C,phi,S12,S1,S2,f]=coherencyc(ZscoredFiberData(:,3),ZscoredFiberData(:,2),params);
% ChunkData.FrequencyDomain.(['Coherence_' OpticalChannelNames{3} '_' OpticalChannelNames{2}])=C;
% ChunkData.FrequencyDomain.(['Phase_' OpticalChannelNames{3} '_' OpticalChannelNames{2}])=phi;
% ChunkData.FrequencyDomain.(['Frequency_' OpticalChannelNames{3} '_' OpticalChannelNames{2}])=f;

% [C,phi,S12,S1,S2,f]=coherencyc(ZscoredFiberData(:,3),ZscoredFiberData(:,1),params);
% ChunkData.FrequencyDomain.(['Coherence_' OpticalChannelNames{3} '_' OpticalChannelNames{1}])=C;
% ChunkData.FrequencyDomain.(['Phase_' OpticalChannelNames{3} '_' OpticalChannelNames{1}])=phi;
% ChunkData.FrequencyDomain.(['Frequency_' OpticalChannelNames{3} '_' OpticalChannelNames{1}])=f;

% [C,phi,S12,S1,S2,f]=coherencyc(ZscoredFiberData(:,2),ZscoredFiberData(:,1),params);
% ChunkData.FrequencyDomain.(['Coherence_' OpticalChannelNames{2} '_' OpticalChannelNames{1}])=C;
% ChunkData.FrequencyDomain.(['Phase_' OpticalChannelNames{2} '_' OpticalChannelNames{1}])=phi;
% ChunkData.FrequencyDomain.(['Frequency_' OpticalChannelNames{2} '_' OpticalChannelNames{1}])=f;

%% Binarize locomotion data
[imp_bin]=velocity_binarize_fiberphotometry(WheelData,ChunkData.Params.DataFs,ChunkData.Params.DataFs,1e-5);

%% Cross correlation analysis
ExpectedLength=(ChunkData.Params.DataFs*ChunkData.Params.DataSeconds)-StartInd;
maxLag=60*ChunkData.Params.DataFs; %Calculate cross correlation with +/- 5s lags

[r_GCaMP_Loco]=xcorr(detrend(ZscoredFiberData((1:ExpectedLength),2),'linear'),detrend(double(imp_bin(1:ExpectedLength)),'linear'),maxLag,'coeff');

[r_CBV_Loco]=xcorr(detrend(ZscoredFiberData((1:ExpectedLength),3),'linear'),detrend(double(imp_bin(1:ExpectedLength)),'linear'),maxLag,'coeff');

[r_CBV_GCaMP,lags]=xcorr(detrend(ZscoredFiberData(:,3),'linear'),detrend(ZscoredFiberData(:,2),'linear'),maxLag,'coeff');
% [r_CBV_GFP,lags]=xcorr(detrend(ZscoredFiberData((1:ExpectedLength),3),'linear'),detrend(ZscoredFiberData((1:ExpectedLength),1),'linear'),maxLag,'coeff');

%% Run triggered averaging
RunDuration=[5 10 15 30 45];
RunLengths={'five_second_events','ten_second_events','fifteen_second_events','thirty_second_events','fortyfive_second_events'};
LeadTime=5;    
[~,~,new_T_run]=motion_cont_fiberphotometry(imp_bin,ChunkData.Params.DataFs,RunDuration(1),LeadTime);
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
        
        ChunkData.WheelData.Locomotion.(RunLabel).RunInds=EventInds;
        ChunkData.WheelData.Locomotion.(RunLabel).RunLengths=(EventInds(2,:)-EventInds(1,:))/ChunkData.Params.DataFs;
        ChunkData.AveragedData.Locomotion.(RunLabel).AvgResp=AverageLocomotionResponse;
        ChunkData.AveragedData.Locomotion.(RunLabel).StdResp=StanDevLocomotion;
        ChunkData.AveragedData.Locomotion.(RunLabel).MedResp=MedianLocomotionResponse;
        ChunkData.AveragedData.Locomotion.(RunLabel).BaseResp=BaseResp;
        ChunkData.AveragedData.Locomotion.(RunLabel).PeakResp=PeakResp;
        ChunkData.AveragedData.Locomotion.(RunLabel).PeakVal=PeakVal;
        ChunkData.AveragedData.Locomotion.(RunLabel).PeakTime=PeakTime;
        ChunkData.LocomotionEvokedData.(RunLabel).OpticalData=LocomotionEvents;
        ChunkData.LocomotionEvokedData.(RunLabel).OpticalDataRaw=LocomotionEventsRaw;
    else
        ChunkData.WheelData.Locomotion.(RunLabel).RunInds=[];
        ChunkData.WheelData.Locomotion.(RunLabel).RunLengths=[];
        ChunkData.AveragedData.Locomotion.(RunLabel).AvgResp=[];
        ChunkData.AveragedData.Locomotion.(RunLabel).StdResp=[];
        ChunkData.AveragedData.Locomotion.(RunLabel).MedResp=[];
        ChunkData.AveragedData.Locomotion.(RunLabel).BaseResp=[];
        ChunkData.AveragedData.Locomotion.(RunLabel).PeakResp=[];
        ChunkData.AveragedData.Locomotion.(RunLabel).PeakVal=[];
        ChunkData.AveragedData.Locomotion.(RunLabel).PeakTime=[];
        ChunkData.LocomotionEvokedData.(RunLabel).OpticalData=[];
    end
    
%     plotTime=((1:length(ChunkData.AveragedData.(RunLabel).AvgResp))-ChunkData.Params.StartPad)/ChunkData.Params.DataFs;
%     figure;plot(plotTime,ChunkData.AveragedData.(RunLabel).AvgResp(:,(2:3))); xlim([-5 RunDuration(durNum)]); legend({'nNosCre GCaMP6s','TRITC Blood Volume'});
%     title('Average locomotion evoked response'); xlabel('Time (sec)'); ylabel('Z score');
    clear LocomotionEvents LocomotionEventsRaw
end

%% ROI defined behavior triggered averaging
roiNames=fieldnames(trackingData.binarizedData);
eventDuration=[2 5 10 15 30 45];
eventLengths={'two_second_events','five_second_events','ten_second_events','fifteen_second_events','thirty_second_events','fortyfive_second_events'};
LeadTime=5;
ChunkData.Params.behavior_camFs=length(trackingData.Tracking.(roiNames{1}).smoothTracking)/ChunkData.Params.DataSeconds;
for roiNum=1:size(roiNames,1)
    if ~isempty(trackingData.binarizedData.(roiNames{roiNum}).behaviorInds)
        EventLengths=(trackingData.binarizedData.(roiNames{roiNum}).behaviorInds(2,:)-trackingData.binarizedData.(roiNames{roiNum}).behaviorInds(1,:))/ChunkData.Params.behavior_camFs;
        behaviorLabel=[roiNames{roiNum} 'EvokedData'];
        for durNum=1:length(eventDuration)
            behaviorEvents=[];
            eventLabel=eventLengths{durNum};
            
            if durNum==length(eventDuration)
                EventFind=EventLengths>=eventDuration(durNum);
                EventInds=round(((trackingData.binarizedData.(roiNames{roiNum}).behaviorInds(:,EventFind)...
                    /ChunkData.Params.behavior_camFs)*ChunkData.Params.DataFs),0); % converts frame indicies to doric Fs indices
            else
                EventFind=EventLengths>=eventDuration(durNum) & EventLengths<eventDuration(durNum+1);
                EventInds=round(((trackingData.binarizedData.(roiNames{roiNum}).behaviorInds(:,EventFind)...
                    /ChunkData.Params.behavior_camFs)*ChunkData.Params.DataFs),0);% converts frame indicies to doric Fs indices
            end
            [row,col]=find(EventInds<=(60*ChunkData.Params.DataFs));
            EventInds(:,col)=[];
            EventInds=EventInds-(60*ChunkData.Params.DataFs);
            for k=1:size(EventInds,2)
                StartInd=EventInds(1,k)-ChunkData.Params.StartPad;
                if durNum==length(eventDuration)
                    EndInd=EventInds(1,k)+(ChunkData.Params.DataFs*60);
                else
                    EndInd=EventInds(1,k)+(ChunkData.Params.DataFs*eventDuration(durNum+1));%ChunkData.Params.FollowPad;
                end
                if StartInd>0
                    if EndInd<length(ZscoredFiberData)+1
                        Baseline=mean(ZscoredFiberData((StartInd:(StartInd+(2*ChunkData.Params.DataFs))),:),1);
                        BaselineRaw=mean(LowPassData((StartInd:(StartInd+(2*ChunkData.Params.DataFs))),:),1);
                        for j=1:size(Baseline,2)
                            behaviorEvents(:,j,k)=ZscoredFiberData((StartInd:EndInd),j)-Baseline(j);
                            behaviorEventsRaw(:,j,k)=LowPassData((StartInd:EndInd),j)-BaselineRaw(j);
                        end
                    end
                end
            end
            
            if ~isempty(behaviorEvents)
                AveragebehaviorResponse=mean(behaviorEvents,3);
                StanDevBehavior=std(behaviorEvents,0,3);
                MedianbehaviorResponse=median(behaviorEvents,3);
                BaseResp=mean(AveragebehaviorResponse((1:(3*ChunkData.Params.DataFs)),:),1);
                PeakVal=max(AveragebehaviorResponse,[],1);
                for chanNum=1:size(AveragebehaviorResponse,2)
                    PeakTime(chanNum)=(find(AveragebehaviorResponse(:,chanNum)==PeakVal(chanNum))-ChunkData.Params.StartPad)/ChunkData.Params.DataFs;
                end
                PeakInds=((ChunkData.Params.StartPad+eventDuration(durNum)*ChunkData.Params.DataFs)-(3*ChunkData.Params.DataFs)):(ChunkData.Params.StartPad+eventDuration(durNum)*ChunkData.Params.DataFs);
                PeakResp=mean(AveragebehaviorResponse(PeakInds,:),1);
                
                ChunkData.WheelData.(roiNames{roiNum}).(eventLabel).RunInds=EventInds;
                ChunkData.WheelData.(roiNames{roiNum}).(eventLabel).RunLengths=(EventInds(2,:)-EventInds(1,:))/ChunkData.Params.DataFs;
                ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).AvgResp=AveragebehaviorResponse;
                ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).StdResp=StanDevBehavior;
                ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).MedResp=MedianbehaviorResponse;
                ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).BaseResp=BaseResp;
                ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).PeakResp=PeakResp;
                ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).PeakVal=PeakVal;
                ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).PeakTime=PeakTime;
                ChunkData.(behaviorLabel).(eventLabel).OpticalData=behaviorEvents;
                ChunkData.(behaviorLabel).(eventLabel).OpticalDataRaw=behaviorEventsRaw;
                ChunkData.(behaviorLabel).roiTracks.rawTrack=trackingData.Tracking.(roiNames{roiNum}).rawTracking;
                ChunkData.(behaviorLabel).roiTracks.smoothTrack=trackingData.Tracking.(roiNames{roiNum}).smoothTracking;
                ChunkData.(behaviorLabel).binData.binTracking=trackingData.binarizedData.(roiNames{roiNum}).binBehavior;
                ChunkData.(behaviorLabel).binData.eventInds=trackingData.binarizedData.(roiNames{roiNum}).behaviorInds;
                ChunkData.(behaviorLabel).binData.binThresh=trackingData.binarizedData.(roiNames{roiNum}).behaviorThresh;
            else
                ChunkData.WheelData.(roiNames{roiNum}).(eventLabel).RunInds=[];
                ChunkData.WheelData.(roiNames{roiNum}).(eventLabel).RunLengths=[];
                ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).AvgResp=[];
                ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).StdResp=[];
                ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).MedResp=[];
                ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).BaseResp=[];
                ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).PeakResp=[];
                ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).PeakVal=[];
                ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).PeakTime=[];
                ChunkData.(behaviorLabel).(eventLabel).OpticalData=[];
                ChunkData.(behaviorLabel).(eventLabel).OpticalDataRaw=[];
                ChunkData.(behaviorLabel).roiTracks.rawTrack=trackingData.Tracking.(roiNames{roiNum}).rawTracking;
                ChunkData.(behaviorLabel).roiTracks.smoothTrack=trackingData.Tracking.(roiNames{roiNum}).smoothTracking;
                ChunkData.(behaviorLabel).binData.binTracking=trackingData.binarizedData.(roiNames{roiNum}).binBehavior;
                ChunkData.(behaviorLabel).binData.eventInds=trackingData.binarizedData.(roiNames{roiNum}).behaviorInds;
                ChunkData.(behaviorLabel).binData.binThresh=trackingData.binarizedData.(roiNames{roiNum}).behaviorThresh;
            end
            
            clear behaviorEvents behaviorEventsRaw
        end
    else
        ChunkData.WheelData.(roiNames{roiNum}).(eventLabel).RunInds=[];
        ChunkData.WheelData.(roiNames{roiNum}).(eventLabel).RunLengths=[];
        ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).AvgResp=[];
        ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).StdResp=[];
        ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).MedResp=[];
        ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).BaseResp=[];
        ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).PeakResp=[];
        ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).PeakVal=[];
        ChunkData.AveragedData.(roiNames{roiNum}).(eventLabel).PeakTime=[];
        ChunkData.(behaviorLabel).(eventLabel).OpticalData=[];
        ChunkData.(behaviorLabel).(eventLabel).OpticalDataRaw=[];
        ChunkData.(behaviorLabel).roiTracks.rawTrack=trackingData.Tracking.(roiNames{roiNum}).rawTracking;
        ChunkData.(behaviorLabel).roiTracks.smoothTrack=trackingData.Tracking.(roiNames{roiNum}).smoothTracking;
        ChunkData.(behaviorLabel).binData.binTracking=trackingData.binarizedData.(roiNames{roiNum}).binBehavior;
        ChunkData.(behaviorLabel).binData.eventInds=trackingData.binarizedData.(roiNames{roiNum}).behaviorInds;
        ChunkData.(behaviorLabel).binData.binThresh=trackingData.binarizedData.(roiNames{roiNum}).behaviorThresh;
    end
end



%% House keeping of raw data
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



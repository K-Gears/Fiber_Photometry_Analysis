function [ChunkData]=ExtractFPdata_Pharmacology_001(filename,OpticalChannelNames,AnalogChannelNames,CorrectionConst,writeTime)
%% Read Me
% Function used to extract imaging data from .csv files generated during
% fiber photometry experiments.
%filename: cell array of 3 filenames covering the three periods of the
%experiment.
%_1: pre injection time period used for establishing baseline values
%_2:post injection wash on period immediately following injection
%_3: post injection steady state period following injection
%OpticalChannelNames:Cell array containing names to be assigned to optical
%channels. {'autofluoresence','EGFP','TRITC'}
%AnalogChannelNames:Cell array containing names to be assigned to analog
%signals {'wheel','left whisker solenoid','right whisker solenoid}

close all
%% Constants
ChunkData.Params.Channels=[OpticalChannelNames,AnalogChannelNames];
directory=cd;
folderBreaks=strfind(directory,'\');
fileBreaks=strfind(filename{1},'_');
ChunkData.Params.animalname=filename{1}(1:9);
ChunkData.Params.AnimalType=filename{1}((fileBreaks(2)+1):fileBreaks(4)-1);
if numel(folderBreaks)==4
    ChunkData.Params.date=directory((folderBreaks(2)+1):(folderBreaks(3)-1));
    ChunkData.Params.Injection_Type=directory((folderBreaks(2)+1):(folderBreaks(3)-1));
    ChunkData.Params.fiber_depth=directory((folderBreaks(3)+1):(folderBreaks(4)-1));    
else
    ChunkData.Params.date=directory((folderBreaks(2)+1):(folderBreaks(3)-1));
    ChunkData.Params.Injection_Type=directory((folderBreaks(3)+1):(folderBreaks(4)-1));
    ChunkData.Params.fiber_depth=directory((folderBreaks(4)+1):(folderBreaks(5)-1));   
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
ChunkData.Params.HbTCorrectionConstant=CorrectionConst;%Factor to multipy TRITC by before adding to GCaMP to remove blood volume effects
ChunkData.Params.FileDuration=[60,20,60];%duration in minutes for each acquisition file
params.Fs=ChunkData.Params.DataFs; %multitaper estimation parameters
params.tapers=[3 5];%multitaper estimation parameters

[z,p,k]=butter(3,ChunkData.Params.Final_Freq/(0.5*ChunkData.Params.DataFs),'bandpass'); %Low pass filter for locomotion data
[sos_final,g_final]=zp2sos(z,p,k);

[z,p,k]=butter(3,10/(0.5*ChunkData.Params.DataFs),'low'); %Low pass filter for locomotion data
[sos_ball,g_ball]=zp2sos(z,p,k);

%% Read .csv files
DataChannels=[2,3,5,7];%405nm,465nm,540nm,wheel analog signal

PreInjectionData=csvread(filename{1},2,0);
PreInjectionData=PreInjectionData(:,DataChannels);
BandpassPre=filtfilt(sos_final,g_final,PreInjectionData);
expectedlength=ChunkData.Params.FileDuration(1)*60*ChunkData.Params.DataFs;
BandpassPre=BandpassPre((1:expectedlength),:);
BandpassPre(:,2)=BandpassPre(:,2)+(CorrectionConst*BandpassPre(:,3)); %Correct for hemodynamic attenuation of GCaMP signal
BandpassPreFreq=BandpassPre;
ClickInds=[];
for chanNum=2:3
    ClickFind=abs(diff(BandpassPre(:,chanNum),1,1));
    ClickInds=[ClickInds;find(ClickFind>=2e-4)+1];
end
ClickInds=unique(ClickInds);
clickCnt=1;
ClickWin=[];
while clickCnt<length(ClickInds)
    if clickCnt==1
        ClickWin(1,clickCnt)=ClickInds(clickCnt);
        ClickWin(2,clickCnt)=ClickInds(clickCnt)+(3*ChunkData.Params.DataFs+1);
        clickCnt=clickCnt+1;
    elseif (ClickInds(clickCnt)-ClickInds((clickCnt-1)))<=(3*ChunkData.Params.DataFs+1)
        ClickInds(clickCnt)=[];
    else
        ClickWin(1,clickCnt)=ClickInds(clickCnt);
        ClickWin(2,clickCnt)=ClickInds(clickCnt)+(3*ChunkData.Params.DataFs+1);
        clickCnt=clickCnt+1;
    end
end
for clickNum=1:size(ClickWin,2)
BandpassPre((ClickWin(1,clickNum):ClickWin(2,clickNum)),(2:3))=NaN;
end


WashOnData=csvread(filename{2},2,0);
WashOnData=WashOnData(:,DataChannels);
BandpassWash=filtfilt(sos_final,g_final,WashOnData);
expectedlength=ChunkData.Params.FileDuration(2)*60*ChunkData.Params.DataFs;
BandpassWash=BandpassWash((1:expectedlength),:);
BandpassWash(:,2)=BandpassWash(:,2)+(CorrectionConst*BandpassWash(:,3)); %Correct for hemodynamic attenuation of GCaMP signal
BandpassWashFreq=BandpassWash;
ClickInds=[];
for chanNum=2:3
    ClickFind=abs(diff(BandpassWash(:,chanNum),1,1));
    ClickInds=[ClickInds;find(ClickFind>=2e-4)+1];
end
ClickInds=unique(ClickInds);
clickCnt=1;
ClickWin=[];
while clickCnt<length(ClickInds)
    if clickCnt==1
        ClickWin(1,clickCnt)=ClickInds(clickCnt);
        ClickWin(2,clickCnt)=ClickInds(clickCnt)+(3*ChunkData.Params.DataFs+1);
        clickCnt=clickCnt+1;
    elseif (ClickInds(clickCnt)-ClickInds((clickCnt-1)))<=(3*ChunkData.Params.DataFs+1)
        ClickInds(clickCnt)=[];
    else
        ClickWin(1,clickCnt)=ClickInds(clickCnt);
        ClickWin(2,clickCnt)=ClickInds(clickCnt)+(3*ChunkData.Params.DataFs+1);
        clickCnt=clickCnt+1;
    end
end
for clickNum=1:size(ClickWin,2)
BandpassWash((ClickWin(1,clickNum):ClickWin(2,clickNum)),(2:3))=NaN;
end


SteadyStateData=csvread(filename{3},2,0);
SteadyStateData=SteadyStateData(:,DataChannels);
if strcmpi(ChunkData.Params.date,'062520') && strcmpi(ChunkData.Params.animalname,'CE_FBR029')
    SteadyStateData(:,3)=SteadyStateData(:,3)*10;
end
BandpassSteady=filtfilt(sos_final,g_final,SteadyStateData);
expectedlength=ChunkData.Params.FileDuration(3)*60*ChunkData.Params.DataFs;
BandpassSteady=BandpassSteady((1:expectedlength),:);
BandpassSteady(:,2)=BandpassSteady(:,2)+(CorrectionConst*BandpassSteady(:,3)); %Correct for hemodynamic attenuation of GCaMP signal
BandpassSteadyFreq=BandpassSteady;
ClickInds=[];
for chanNum=2:3
    ClickFind=abs(diff(BandpassSteady(:,chanNum),1,1));
    ClickInds=[ClickInds;find(ClickFind>=2e-4)+1];
end
ClickInds=unique(ClickInds);
clickCnt=1;
ClickWin=[];
while clickCnt<length(ClickInds)
    if clickCnt==1
        ClickWin(1,clickCnt)=ClickInds(clickCnt);
        ClickWin(2,clickCnt)=ClickInds(clickCnt)+(3*ChunkData.Params.DataFs+1);
        clickCnt=clickCnt+1;
    elseif (ClickInds(clickCnt)-ClickInds((clickCnt-1)))<=(3*ChunkData.Params.DataFs+1)
        ClickInds(clickCnt)=[];
    else
        ClickWin(1,clickCnt)=ClickInds(clickCnt);
        ClickWin(2,clickCnt)=ClickInds(clickCnt)+(3*ChunkData.Params.DataFs+1);
        clickCnt=clickCnt+1;
    end
end
for clickNum=1:size(ClickWin,2)
BandpassSteady((ClickWin(1,clickNum):ClickWin(2,clickNum)),(2:3))=NaN;
end

%% Z-score optical preinjection data
OffsetStart=120*ChunkData.Params.DataFs; %Removes filtering artifact
AvgData=nanmean(BandpassPre((OffsetStart:end),(1:3)),1);
StdData=std(BandpassPre((OffsetStart:end),(1:3)),0,1,'omitnan');
AvgMatrix=repmat(AvgData,length(BandpassPre((OffsetStart:end),:)),1);
StdMatrix=repmat(StdData,length(BandpassPre((OffsetStart:end),:)),1);
ZscoredPre=(BandpassPre((OffsetStart:end),(1:3))-AvgMatrix)./StdMatrix;
ZscoredPreFreq=(BandpassPreFreq((OffsetStart:end),(1:3))-AvgMatrix)./StdMatrix;

ZscoredSteady=(BandpassSteady((OffsetStart:end),(1:3))-AvgMatrix)./StdMatrix;
ZscoredSteadyFreq=(BandpassSteadyFreq((OffsetStart:end),(1:3))-AvgMatrix)./StdMatrix;

AvgMatrix=repmat(AvgData,length(BandpassWash((OffsetStart:end),:)),1);
StdMatrix=repmat(StdData,length(BandpassWash((OffsetStart:end),:)),1);
ZscoredWash=(BandpassWash((OffsetStart:end),(1:3))-AvgMatrix)./StdMatrix;
ZscoredWashFreq=(BandpassWashFreq((OffsetStart:end),(1:3))-AvgMatrix)./StdMatrix;

%% Concatenate Files
trialPad((1:OffsetStart),(1:length(DataChannels)-1))=NaN;

dTa=caldiff(writeTime(1:2),'Time');
padDura=seconds(split(dTa,'time'))-(ChunkData.Params.FileDuration(2)*60);
padSampa=padDura*ChunkData.Params.DataFs;
padMATa((1:padSampa),(1:length(DataChannels)-1))=NaN;

dTb=caldiff(writeTime(2:3),'Time');
padDurb=seconds(split(dTb,'time'))-(ChunkData.Params.FileDuration(3)*60);
padSampb=padDurb*ChunkData.Params.DataFs;
padMATb((1:padSampb),(1:length(DataChannels)-1))=NaN;

catData=[];
catData=[ZscoredPre;padMATa];
catData=[catData;trialPad];
catData=[catData;ZscoredWash];
catData=[catData;padMATb];
catData=[catData;trialPad];
allData=[catData;ZscoredSteady];
catData=[];

ChunkData.TrialData.ZScored.concatenatedTrial=allData;
ChunkData.TrialData.ZScored.Preinjection=ZscoredPre;
ChunkData.TrialData.ZScored.WashOn=ZscoredWash;
ChunkData.TrialData.ZScored.SteadyState=ZscoredSteady;
ChunkData.TrialData.RawSignal.preInjection=BandpassPre;
ChunkData.TrialData.RawSignal.washOn=BandpassWash;
ChunkData.TrialData.RawSignal.steadyState=BandpassSteady;


%% Cross Spectral Frequency Coherence
params.fpass=[0.1 1];
params.tapers=[19 37];

[C,phi,S12,S1,S2,f]=coherencyc(ZscoredPreFreq(:,3),ZscoredPreFreq(:,2),params);
ChunkData.FrequencyDomain.PreInjection.(['Coherence_' OpticalChannelNames{3} '_' OpticalChannelNames{2}])=C;
ChunkData.FrequencyDomain.PreInjection.(['Phase_' OpticalChannelNames{3} '_' OpticalChannelNames{2}])=phi;
ChunkData.FrequencyDomain.PreInjection.(['Frequency_' OpticalChannelNames{3} '_' OpticalChannelNames{2}])=f;

[C,phi,S12,S1,S2,f]=coherencyc(ZscoredSteadyFreq(:,3),ZscoredSteadyFreq(:,2),params);
ChunkData.FrequencyDomain.SteadyState.(['Coherence_' OpticalChannelNames{3} '_' OpticalChannelNames{2}])=C;
ChunkData.FrequencyDomain.SteadyState.(['Phase_' OpticalChannelNames{3} '_' OpticalChannelNames{2}])=phi;
ChunkData.FrequencyDomain.SteadyState.(['Frequency_' OpticalChannelNames{3} '_' OpticalChannelNames{2}])=f;

%% Cross correlation analysis
maxLag=60*ChunkData.Params.DataFs; %Calculate cross correlation with +/- 5s lag

[r_CBV_GCaMP]=xcorr(detrend(ZscoredPreFreq(:,3),'linear'),detrend(ZscoredPreFreq(:,2),'linear'),maxLag,'coeff');
ChunkData.Xcorr.PreInjection.CBV_GCaMP=r_CBV_GCaMP;

[r_CBV_GCaMP,lags]=xcorr(detrend(ZscoredSteadyFreq(:,3),'linear'),detrend(ZscoredSteadyFreq(:,2),'linear'),maxLag,'coeff');
ChunkData.Xcorr.SteadyState.CBV_GCaMP=r_CBV_GCaMP;
ChunkData.Xcorr.Lags=lags/ChunkData.Params.DataFs;

%% Format data for run triggered averaging
StateLabel={'PreInjection','SteadyState'};
WheelData(:,1)=filtfilt(sos_ball,g_ball,PreInjectionData((1:expectedlength),4));
WheelData(:,2)=filtfilt(sos_ball,g_ball,SteadyStateData((1:expectedlength),4));
ZscoreData(:,:,1)=ZscoredPre;
ZscoreData(:,:,2)=ZscoredSteady;
for stateNum=1:2
%% Binarize locomotion data
[imp_bin]=velocity_binarize(WheelData((OffsetStart:end),stateNum),ChunkData.Params.DataFs,ChunkData.Params.DataFs,1e-5);

%% Run triggered averaging
RunDuration=[5 10 15 30 45];
RunLengths={'five_second_events','ten_second_events','fifteen_second_events','thirty_second_events','fortyfive_second_events'};
LeadTime=5;    
[~,~,new_T_run]=motion_cont_2(imp_bin,ChunkData.Params.DataFs,RunDuration(1),LeadTime);
EventLengths=(new_T_run(2,:)-new_T_run(1,:))/ChunkData.Params.DataFs;
for eventNum=1:size(new_T_run,2)
    StartInd=new_T_run(1,eventNum)-ChunkData.Params.StartPad;
    EndInd=new_T_run(1,eventNum)+ChunkData.Params.StartPad;
    Baseline=mean(ZscoreData((StartInd:(StartInd+(2*ChunkData.Params.DataFs))),:,stateNum),1);
    for j=1:size(Baseline,2)
        LocomotionEvents(:,j,eventNum)=ZscoreData((StartInd:EndInd),j,stateNum)-Baseline(j);
    end
end
if ~isempty(LocomotionEvents)
    AverageLocomotionResponse=nanmean(LocomotionEvents,3);
    StanDevLocomotion=std(LocomotionEvents,0,3,'omitnan');
    MedianLocomotionResponse=median(LocomotionEvents,3);
    BaseResp=nanmean(AverageLocomotionResponse((1:(3*ChunkData.Params.DataFs)),:),1);
    PeakInds=((ChunkData.Params.StartPad+RunDuration(1)*ChunkData.Params.DataFs)-(3*ChunkData.Params.DataFs)):(ChunkData.Params.StartPad+RunDuration(1)*ChunkData.Params.DataFs);
    PeakResp=nanmean(AverageLocomotionResponse(PeakInds,:),1);
    
%     ChunkData.WheelData.(StateLabel{stateNum}).All.RunInds=EventInds;
    ChunkData.WheelData.(StateLabel{stateNum}).All.RunLengths=(new_T_run(2,:)-new_T_run(1,:))/ChunkData.Params.DataFs;
    ChunkData.AveragedData.(StateLabel{stateNum}).All.AvgResp=AverageLocomotionResponse;
    ChunkData.AveragedData.(StateLabel{stateNum}).All.StdResp=StanDevLocomotion;
    ChunkData.AveragedData.(StateLabel{stateNum}).All.MedResp=MedianLocomotionResponse;
    ChunkData.AveragedData.(StateLabel{stateNum}).All.BaseResp=BaseResp;
    ChunkData.AveragedData.(StateLabel{stateNum}).All.PeakResp=PeakResp;
    ChunkData.AveragedData.(StateLabel{stateNum}).All.RespVol=sum(AverageLocomotionResponse((ChunkData.Params.StartPad:end),:),1)/length(AverageLocomotionResponse((ChunkData.Params.StartPad:end),1));
    ChunkData.LocomotionEvokedData.(StateLabel{stateNum}).All.OpticalData=LocomotionEvents;
else
    ChunkData.WheelData.(StateLabel{stateNum}).All.RunInds=[];
    ChunkData.WheelData.(StateLabel{stateNum}).All.RunLengths=[];
    ChunkData.AveragedData.(StateLabel{stateNum}).All.AvgResp=[];
    ChunkData.AveragedData.(StateLabel{stateNum}).All.StdResp=[];
    ChunkData.AveragedData.(StateLabel{stateNum}).All.MedResp=[];
    ChunkData.AveragedData.(StateLabel{stateNum}).All.BaseResp=[];
    ChunkData.AveragedData.(StateLabel{stateNum}).All.PeakResp=[];
    ChunkData.AveragedData.(StateLabel{stateNum}).All.RespVol=[];
    ChunkData.LocomotionEvokedData.(StateLabel{stateNum}).All.OpticalData=[];
end
clear LocomotionEvents
    
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
            if EndInd<length(ZscoreData(:,:,stateNum))+1
                Baseline=mean(ZscoreData((StartInd:(StartInd+(2*ChunkData.Params.DataFs))),:,stateNum),1);
                for j=1:size(Baseline,2)
                    LocomotionEvents(:,j,k)=ZscoreData((StartInd:EndInd),j,stateNum)-Baseline(j);
                end
            end
        end
    end
    
    if ~isempty(LocomotionEvents)
        AverageLocomotionResponse=nanmean(LocomotionEvents,3);
        StanDevLocomotion=std(LocomotionEvents,0,3,'omitnan');
        MedianLocomotionResponse=median(LocomotionEvents,3);
        BaseResp=nanmean(AverageLocomotionResponse((1:(3*ChunkData.Params.DataFs)),:),1);
        PeakInds=((ChunkData.Params.StartPad+RunDuration(durNum)*ChunkData.Params.DataFs)-(3*ChunkData.Params.DataFs)):(ChunkData.Params.StartPad+RunDuration(durNum)*ChunkData.Params.DataFs);
        PeakResp=nanmean(AverageLocomotionResponse(PeakInds,:),1);
        
        ChunkData.WheelData.(StateLabel{stateNum}).(RunLabel).RunInds=EventInds;
        ChunkData.WheelData.(StateLabel{stateNum}).(RunLabel).RunLengths=(EventInds(2,:)-EventInds(1,:))/ChunkData.Params.DataFs;
        ChunkData.AveragedData.(StateLabel{stateNum}).(RunLabel).AvgResp=AverageLocomotionResponse;
        ChunkData.AveragedData.(StateLabel{stateNum}).(RunLabel).StdResp=StanDevLocomotion;
        ChunkData.AveragedData.(StateLabel{stateNum}).(RunLabel).MedResp=MedianLocomotionResponse;
        ChunkData.AveragedData.(StateLabel{stateNum}).(RunLabel).BaseResp=BaseResp;
        ChunkData.AveragedData.(StateLabel{stateNum}).(RunLabel).PeakResp=PeakResp;
        ChunkData.AveragedData.(StateLabel{stateNum}).(RunLabel).RespVol=sum(AverageLocomotionResponse((ChunkData.Params.StartPad:end),:),1)/length(AverageLocomotionResponse((ChunkData.Params.StartPad:end),1));
        ChunkData.LocomotionEvokedData.(StateLabel{stateNum}).(RunLabel).OpticalData=LocomotionEvents;
    else
        ChunkData.WheelData.(StateLabel{stateNum}).(RunLabel).RunInds=[];
        ChunkData.WheelData.(StateLabel{stateNum}).(RunLabel).RunLengths=[];
        ChunkData.AveragedData.(StateLabel{stateNum}).(RunLabel).AvgResp=[];
        ChunkData.AveragedData.(StateLabel{stateNum}).(RunLabel).StdResp=[];
        ChunkData.AveragedData.(StateLabel{stateNum}).(RunLabel).MedResp=[];
        ChunkData.AveragedData.(StateLabel{stateNum}).(RunLabel).BaseResp=[];
        ChunkData.AveragedData.(StateLabel{stateNum}).(RunLabel).PeakResp=[];
        ChunkData.AveragedData.(StateLabel{stateNum}).(RunLabel).RespVol=[];
        ChunkData.LocomotionEvokedData.(StateLabel{stateNum}).(RunLabel).OpticalData=[];
    end
    clear LocomotionEvents
end
end
plotTime=((1:length(ChunkData.AveragedData.PreInjection.All.AvgResp))-ChunkData.Params.StartPad)/ChunkData.Params.DataFs;
figure;hold on;

plot(plotTime,ChunkData.AveragedData.PreInjection.All.AvgResp(:,2),'Color','b'); %locomotion evoked average
plot(plotTime,(ChunkData.AveragedData.PreInjection.All.AvgResp(:,2)+ChunkData.AveragedData.PreInjection.All.StdResp(:,2)),'LineStyle','--','Color','b','HandleVisibility','off');
plot(plotTime,(ChunkData.AveragedData.PreInjection.All.AvgResp(:,2)-ChunkData.AveragedData.PreInjection.All.StdResp(:,2)),'LineStyle','--','Color','b','HandleVisibility','off');

plot(plotTime,ChunkData.AveragedData.SteadyState.All.AvgResp(:,2),'Color','g'); %locomotion evoked average
plot(plotTime,(ChunkData.AveragedData.SteadyState.All.AvgResp(:,2)+ChunkData.AveragedData.SteadyState.All.StdResp(:,2)),'LineStyle','--','Color','g','HandleVisibility','off');
plot(plotTime,(ChunkData.AveragedData.SteadyState.All.AvgResp(:,2)-ChunkData.AveragedData.SteadyState.All.StdResp(:,2)),'LineStyle','--','Color','g','HandleVisibility','off');

title('Locomotion evoked GCaMP6s intensity change')
ylabel('Z units');
xlabel('Time (s)');
legend({'Pre-injection',ChunkData.Params.Injection_Type});

figure;hold on;

plot(plotTime,ChunkData.AveragedData.PreInjection.All.AvgResp(:,3),'Color','b'); %locomotion evoked average
plot(plotTime,(ChunkData.AveragedData.PreInjection.All.AvgResp(:,3)+ChunkData.AveragedData.PreInjection.All.StdResp(:,3)),'LineStyle','--','Color','b','HandleVisibility','off');
plot(plotTime,(ChunkData.AveragedData.PreInjection.All.AvgResp(:,3)-ChunkData.AveragedData.PreInjection.All.StdResp(:,3)),'LineStyle','--','Color','b','HandleVisibility','off');

plot(plotTime,ChunkData.AveragedData.SteadyState.All.AvgResp(:,3),'Color','g'); %locomotion evoked average
plot(plotTime,(ChunkData.AveragedData.SteadyState.All.AvgResp(:,3)+ChunkData.AveragedData.SteadyState.All.StdResp(:,3)),'LineStyle','--','Color','g','HandleVisibility','off');
plot(plotTime,(ChunkData.AveragedData.SteadyState.All.AvgResp(:,3)-ChunkData.AveragedData.SteadyState.All.StdResp(:,3)),'LineStyle','--','Color','g','HandleVisibility','off');

title('Locomotion evoked TRITC intensity change')
ylabel('Z units');
xlabel('Time (s)');
legend({'Pre-injection',ChunkData.Params.Injection_Type});

save([ChunkData.Params.animalname '_' ChunkData.Params.date '_' ChunkData.Params.Injection_Type '_FiberPhotometry.mat'],'ChunkData','-v7.3');

end


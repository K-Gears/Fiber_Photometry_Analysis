function [binBehavior,new_eventInds,movementThresh]=binarize_BehaviorCam(movementData,restDur,eventDur,roiName,movementThresh,varargin)
%This function converts ROI calculated movement data in to binarized
%behaving(1) quiescent(0) bins based on a user defined threshold and
%returns the indicies of events using user defined event length duration
%and preceding periods of rest
% INPUTS:
%movementData:vector of ROI defined movement
%movementThresh: threshold to binarize movement data
%restDur: time in samples to precede and movement event for chunking data
%eventDur:time in samples for an event to be included in chunking data
%
%Outputs:
%binData: logical vector same length as movementData declaring movement above(1) or below(0) movementThresh
%eventInds: 2x number of event vector with top row being event start
%indicies and bottom row being stop indicies.


%% Define threshold if not an input
if ~exist('movementThresh','var')
    threshok='n';
    while strcmpi(threshok,'n')
        close all
        Fs=30;
        plotTime=(1:length(movementData))/Fs;
        figure; hold on;plot(plotTime,movementData,'k');
        xlim([0 7200]); xlabel('Time (s)');
        title(['Select acceleration threshold for binarizing movement of ' roiName '.']);
        movementThresh=input('Define threshold for binarizing movement');
        plotThresh(1:length(plotTime))=movementThresh;
        plot(plotTime,plotThresh,'--r');
        threshok=input('Is threshold ok? (y/n)','s');
    end
end

%% Binarize behavior
Fs=30;
binBehavior=movementData>=movementThresh;

%% Fuse events of rest and behavior below specified durations
[eventInds,restInds]=getInds(binBehavior,movementData);
new_binBehavior=binBehavior;
fillRest=1*Fs;
for restNum=1:size(restInds,2)
    if (restInds(2,restNum)-restInds(1,restNum))<fillRest
        new_binBehavior(restInds(1,restNum):restInds(2,restNum))=1;
    end
end
[eventInds,restInds]=getInds(new_binBehavior,movementData);
final_binBehavior=new_binBehavior;
fillTwitch=0.5*Fs;
for eventNum=1:size(eventInds,2)
    if (eventInds(2,eventNum)-eventInds(1,eventNum))<fillTwitch
        final_binBehavior(eventInds(1,eventNum):eventInds(2,eventNum))=0;
    end
end
[eventInds,restInds]=getInds(final_binBehavior,movementData);

%% Get event inds to be used for behavior chunking
new_eventInds=[];
eventCount=1;
if restInds(1,1)<eventInds(1,1)
    for eventNum=1:size(eventInds,2)
        if (restInds(2,eventNum)-restInds(1,eventNum))>=restDur
            if (eventInds(2,eventNum)-eventInds(1,eventNum))>=eventDur
                new_eventInds(:,eventCount)=eventInds(:,eventNum);
                eventCount=eventCount+1;
            end
        end
    end
else
    for eventNum=2:size(eventInds,2)
        if (restInds(2,(eventNum-1))-restInds(1,(eventNum-1)))>=restDur
            if (eventInds(2,eventNum)-eventInds(1,eventNum))>=eventDur
                new_eventInds(:,eventCount)=eventInds(:,eventNum);
                eventCount=eventCount+1;
            end
        end
    end
end
        
                
    
end

%% Get behavior Indicies
function [eventInds,restInds]=getInds(binBehavior,movementData)
binDiff=diff(binBehavior);

StartInds=find(binDiff==1)+1;
StopInds=find(binDiff==-1);

if StartInds(1)<=StopInds(1)
    if length(StartInds)==length(StopInds)
        eventInds=[StartInds;StopInds];
        restInds(:,1)=[1;(StartInds(1)-1)];
        restInds(:,(2:length(StopInds)))=[(StopInds(1:(end-1))+1);StartInds(2:end)-1];
        restInds(:,(length(StopInds)+1))=[(StopInds(end)+1);length(movementData)];
    else
        eventInds=[StartInds(1:(end-1));StopInds];
        eventInds(:,length(StartInds))=[StartInds(end);length(movementData)];
        restInds(:,1)=[1;(StartInds(1)-1)];
        restInds(:,(2:length(StartInds)))=[(StopInds+1);(StartInds(2:end)-1)];
    end
else
    if length(StartInds)==length(StopInds)
        eventInds(:,1)=[1,StopInds(1)];
        eventInds(:,(2:(length(StopInds))))=[StartInds(1:(end-1));StopInds(2:end)];
        eventInds(:,(length(StopInds)+1))=[StartInds(end);length(movementData)];
        restInds=[(StopInds+1);(StartInds-1)];
    else
        eventInds(:,1)=[1,StopInds(1)];
        eventInds(:,(2:(length(StopInds))))=[StartInds;StopInds(2:end)];
        restInds(:,(1:length(StartInds)))=[(StopInds(1:(end-1))+1);(StartInds-1)];
        restInds(:,length(StopInds))=[(StopInds(end)+1);length(movementData)];
    end
end
end

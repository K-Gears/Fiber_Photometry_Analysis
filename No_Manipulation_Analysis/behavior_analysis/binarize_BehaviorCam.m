function [binData,eventInds,movementThresh]=binarize_BehaviorCam(movementData,movementThresh,restDur,eventDur,varargin)
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
if ~exist(movementThresh,'var')
    Fs=30;
    plotTime=(1:length(movementData))/Fs;
    figure;plot(plotTime,movementData);
    xlim([0 7200]); xlabel('Time (s)');
    movementThresh=input('Define threshold for binarizing movement');
end

%% Binarize behavior
binBehavior=movementData>=movementThresh;

%% Get behavior Inds
binDiff=diff(binBehavior);

StartInds=find(binDiff)+1;
StopInds=find(binDiff);

if StartInds(1)<StopInds(1)
    if length(StartInds)==length(StopInds)
        eventInds=[StartInds;StopInds];
        RestInds(:,1)=[1;StartInds(1)];
        RestInds(:,(2:length(StopInds)))=[(StopInds(1:(end-1))+1);StartInds(2:end)-1];
        RestInds(:,(length(StopInds)+1))=[StopInds(end);length(movementData)];
    else
        eventInds=[StartInds(1:(end-1));StopInds];
        eventInds(:,length(StartInds))=[StartInds(end);length(movementData)];
    end
else
    if length(StartInds)==length(StopInds)
        eventInds(:,1)=[1,StopInds(1)];
        eventInds(:,(2:(length(StopInds))))=[StartInds(1:(end-1));StopInds(2:end)];
        eventInds(:,(length(StopInds)+1))=[StartInds(end);length(movementData)];
    else
        eventInds(:,1)=[1,StopInds(1)];
        eventInds(:,(2:(length(StopInds))))=[StartInds;StopInds(2:end)];
    end

end
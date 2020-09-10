function [BinEdges,NormHist,theFit]=Correction_Histogram(ChunkData)
close all
% BinEdges=(-1:0.002:1);
BinEdges=(-1:0.02:1);
figure(100);
DataHist=histogram2(ChunkData.RawFiberData.RescaledData(:,3),ChunkData.RawFiberData.RescaledData(:,2),BinEdges,BinEdges);
% DataHist=histogram2(ChunkData.RawFiberData.LowPassData(:,3),ChunkData.RawFiberData.LowPassData(:,2),BinEdges,BinEdges);
BinCounts=DataHist.BinCounts;
ColCounts=sum(BinCounts,2);
NormMat=repmat(ColCounts,1,size(BinCounts,2));
NormHist=(BinCounts./NormMat)*100;
KeepBins=ColCounts>=2*ChunkData.Params.DataFs;

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

figure(105); hold on;
imagesc(BinEdges,BinEdges,NormHist');
caxis([0 10]);
h=colorbar('eastoutside');
axis xy;
ylabel('Observed GFP');
xlabel('Observed CBV');
title('Normalized histogram of CBV vs GFP');
plot(BinEdges,PredictedPeaks,'r','LineWidth',2);
% xlim([-0.15 0.25]);
% ylim([-0.1 0.1]);
xlim([0 1]);
ylim([0 1]);

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
theFit=coeffVals(1).*BinEdges+coeffVals(2);
plot(BinEdges,theFit,'c','LineWidth',2);
legend({'Peak count at observed \DeltaCBV', 'Linear fit, peak column values'});
set(get(h,'label'),'string','Percentage of all column counts in bin (%)');
savefig('CE_FBR001_082119_NormalizedHistogramFit');
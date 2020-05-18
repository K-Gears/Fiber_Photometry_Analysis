function plotFig_CE_NOS_Xcorr_005(FigureData)

dataTypes=fieldnames(FigureData.Xcorr.CrossCorr);

MarkerSize=100;
xOffset=[-0.22 0 0.22];

theDate=datetime;
theDate=char(theDate);
theDate=strrep(theDate,'-','');
theDate=strrep(theDate,':','');
theDate=strrep(theDate,' ','_');

for dataNum=1:size(dataTypes,1)

xVals=[];
xVals(1:size(FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).PeakTime,2))=dataNum;
lgndTxt{dataNum}=strrep(dataTypes{dataNum},'_',' ');     

%% Cross-correlation CBV-GCaMP

figure(11); hold on
plot(FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).Lags,FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).CorrCoeff,'Color',FigureData.ColorMaps.cmap((dataNum+1),:),'LineStyle','-','Marker','none');

if dataNum==length(dataTypes)
    xlabel('Lags (sec)');
    ylabel('Corr Coeff');
    title('Cross correlation blood volume  vs GCaMP');
    legend(lgndTxt);
    %     savefig(gcf,['Population_Xcorr_CBV_GFP_' theDate]);
    %     saveas(gcf,['Population_Xcorr_CBV_GFP_' theDate],'epsc');
end

figure(20); hold on
plot(FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).Lags,FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).CorrCoeff,'Color',FigureData.ColorMaps.cmap((dataNum+1),:),'LineStyle','-','Marker','none');
xlim([-5 5]);
if dataNum==length(dataTypes)
    xlabel('Lags (sec)');
    ylabel('Corr Coeff');
    title('Cross correlation blood volume  vs GCaMP');
    legend(lgndTxt);
    %     savefig(gcf,['Population_Xcorr_CBV_GFP_' theDate]);
    %     saveas(gcf,['Population_Xcorr_CBV_GFP_' theDate],'epsc');
end

figure(12); hold on
bar(xVals(1),FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).AvgTime,'FaceColor',FigureData.ColorMaps.cmap((dataNum+1),:),'FaceAlpha',0.3,'EdgeColor','k','LineWidth',1);
scatter(xVals,FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).PeakTime,MarkerSize,'o','filled','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',FigureData.ColorMaps.cmap((dataNum+1),:));

if dataNum==length(dataTypes)
    xlim([0 (size(dataTypes,1)+1)]); xticks(1:1:size(dataTypes,1)); xticklabels(lgndTxt);
    ylabel('Seconds');
    title('Peak correlation lag ');
    %     savefig(gcf,['Population_Peak_Lag_Xcorr_' theDate]);
    %     saveas(gcf,['Population_Peak_Lag_Xcorr_' theDate],'epsc');
end

figure(13); hold on
bar(xVals(1),FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).AvgVal,'FaceColor',FigureData.ColorMaps.cmap((dataNum+1),:),'FaceAlpha',0.3,'EdgeColor','k','LineWidth',1);
scatter(xVals,FigureData.Xcorr.CrossCorr.(dataTypes{dataNum}).PeakVal,MarkerSize,'o','filled','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',FigureData.ColorMaps.cmap((dataNum+1),:));

if dataNum==size(dataTypes,1)
    xlim([0 (size(dataTypes,1)+1)]); xticks(1:1:size(dataTypes,1)); xticklabels(lgndTxt);
    ylim([0 1]);
    ylabel('Corr Coeff');
    title('Maximum Correlation Coefficient');
    %     savefig(gcf,['Population_Peak_Coeff_Xcorr_' theDate]);
    %     saveas(gcf,['Population_Peak_Coeff_Xcorr_' theDate],'epsc');
end

clear PeakTime PeakVal xvals;

%% Coherence CBV-GCaMP

figure(14);set(gca,'XScale','log'); hold on;
semilogx(FigureData.Xcorr.Coherence.(dataTypes{dataNum}).Freqs,FigureData.Xcorr.Coherence.(dataTypes{dataNum}).Coherence,'Color',FigureData.ColorMaps.cmap((dataNum+1),:),'LineStyle','-','Marker','none');

if dataNum==length(dataTypes)
    xlabel('Frequency(Hz)');
    ylabel('Coherence squared');
    title('Coherence blood volume  vs GCaMP');
    legend(lgndTxt);
    %     savefig(gcf,['Population_Coherence_LogFreq_CBV_GFP_' theDate]);
    %     saveas(gcf,['Population_Coherence_LogFreq_CBV_GFP_' theDate],'epsc');
end

end    
clear dataTypes
    %% Signal to Noise
    
    for barNum=1:size(FigureData.Xcorr.SNR.CBVAvg,1)
        xVals=barNum;
        xPos=barNum;
        
        figure(15);hold on;
        bar(xVals,FigureData.Xcorr.SNR.CBVAvg(barNum,:),0.4,'FaceColor',FigureData.ColorMaps.cmap((barNum+1),:),'FaceAlpha',0.3,'EdgeColor','k','LineWidth',0.5);
        
        figure(16);hold on; set(gca,'YScale','log');
        bar(xPos,FigureData.Xcorr.SNR.CBVSTD(barNum,:),0.4,'FaceColor',FigureData.ColorMaps.cmap((barNum+1),:),'FaceAlpha',0.3,'EdgeColor','k','LineWidth',0.5);
        
        figure(17);hold on; set(gca,'YScale','log');
        bar(xPos,FigureData.Xcorr.SNR.GCaMPSTD(barNum,:),0.4,'FaceColor',FigureData.ColorMaps.cmap((barNum+1),:),'FaceAlpha',0.3,'EdgeColor','k','LineWidth',0.5);        
        
        figure(18);hold on;
        bar(xVals,FigureData.Xcorr.SNR.GCaMPAvg(barNum,:),0.4,'FaceColor',FigureData.ColorMaps.cmap((barNum+1),:),'FaceAlpha',0.3,'EdgeColor','k','LineWidth',0.5);
        
    end
    
    fieldTypes=fieldnames(FigureData.Xcorr.SNR);
    d=1;
    for k=1:size(fieldTypes,1)
        findstructs=isstruct(FigureData.Xcorr.SNR.(fieldTypes{k}));
        if findstructs==1
            dataTypes{d,1}=fieldTypes{k};
            d=d+1;
        end
    end
    
    for dataNum=1:size(dataTypes,1)
        titleTxt=strrep(dataTypes{dataNum},'_',' ');
        lgndNames{dataNum}=strrep(dataTypes{dataNum},'_',' ');

        figure(19);hold on;
        plot((0:0.001:7),FigureData.Xcorr.SNR.(dataTypes{dataNum}).GenotypeRawNormPDF,'LineWidth',1,'Color',FigureData.ColorMaps.cmap((dataNum+1),:));

        if dataNum==length(dataTypes)
            xlim([0 1.5]);
            xlabel('Average GCaMP Intensity (V)');
            ylabel('Normalized counts (a.u.)');
            title('Normal distribution of GCaMP intensity values');
            %             savefig(gcf,['Genotype_GCaMP_Raw_Intensity_Distribution_' theDate]);
            %             saveas(gcf,['Genotype_GCaMP_Raw_Intensity_Distribution_' theDate],'epsc');
        end
        
        
        xVals=[];
        xPos=[];
        xVals(1:size(FigureData.Xcorr.SNR.(dataTypes{dataNum}).RawCBVAvg,2))=dataNum;
        xPos(1:size(FigureData.Xcorr.SNR.(dataTypes{dataNum}).RawCBVAvg,2))=dataNum;
        
        figure(15);
        scatter(xVals,FigureData.Xcorr.SNR.(dataTypes{dataNum}).RawCBVAvg,MarkerSize,'o','filled','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',FigureData.ColorMaps.cmap((dataNum+1),:));
        
        if dataNum==size(dataTypes,1)
            xlim([0 (size(dataTypes,1)+1)]); xticks(1:1:size(dataTypes,1)); xticklabels(lgndNames);
            ylabel('Volts');
            title('Blood volume average raw intensity');
%             savefig(['Population_CBV_Avg_Raw_Intensity_' theDate]);
%             saveas(gcf,['Population_CBV_Avg_Raw_Intensity_' theDate],'epsc');
        end
        
        figure(16); set(gca,'YScale','log');
        scatter(xPos,FigureData.Xcorr.SNR.(dataTypes{dataNum}).RawCBVSTD,MarkerSize,'o','filled','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',FigureData.ColorMaps.cmap((dataNum+1),:));
        
        if dataNum==size(dataTypes,1)
            xlim([0 (size(dataTypes,1)+1)]); xticks(1:1:size(dataTypes,1)); xticklabels(lgndNames);
            ylabel('Volts');
            title('Blood volume standard deviation of raw intensity');
%             savefig(['Population_CBV_Raw_Intensity_StandardDev_' theDate]);
%             saveas(gcf,['Population_CBV_Raw_Intensity_StandardDev_' theDate],'epsc');
        end
        
        figure(18);scatter(xVals,FigureData.Xcorr.SNR.(dataTypes{dataNum}).RawGCaMPAvg,MarkerSize,'o','filled','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',FigureData.ColorMaps.cmap((dataNum+1),:));
        if dataNum==size(dataTypes,1)
            xlim([0 (size(dataTypes,1)+1)]); xticks(1:1:size(dataTypes,1)); xticklabels(lgndNames);
            ylabel('Volts');
            title('GCaMP6s average raw intensity');
%             savefig(['Population_GCaMP_Avg_Raw_Intensity_' theDate]);
%             saveas(gcf,['Population_GCaMP_Avg_Raw_Intensity_' theDate],'epsc');
        end
           
        figure(17);set(gca,'YScale','log');
        scatter(xPos,FigureData.Xcorr.SNR.(dataTypes{dataNum}).RawGCaMPSTD,MarkerSize,'o','filled','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',FigureData.ColorMaps.cmap((dataNum+1),:));
        
        if dataNum==size(dataTypes,1)
            xlim([0 (size(dataTypes,1)+1)]); xticks(1:1:size(dataTypes,1)); xticklabels(lgndNames);
            ylabel('Z-units');
            title('GCaMP6s standard deviation of raw intensity');
%             savefig(['Population_GCaMP_Raw_Intensity_StandardDev_' theDate]);
%             saveas(gcf,['Population_GCaMP_Raw_Intensity_StandardDev_' theDate],'epsc');
        end
        
    end
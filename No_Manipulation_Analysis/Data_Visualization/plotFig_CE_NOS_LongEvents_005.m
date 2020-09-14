function plotFig_CE_NOS_LongEvents_005(FigureData)

dataTypes=fieldnames(FigureData.ExtendedEvents.TenSecEventAvg);

MarkerSize=100;
theDate=datetime;
theDate=char(theDate);
theDate=strrep(theDate,'-','');
theDate=strrep(theDate,':','');
theDate=strrep(theDate,' ','_');
xOffset=[-0.22 0 0.22];
thegreens={'#1A9850'; '#66BD63'; '#006837'};
thereds={'#B30000'; '#D7301F'; '#7F0000'};


figure(31); hold on
barHandle=bar(FigureData.ExtendedEvents.BarPlots.CBVAvg','FaceColor','flat','FaceAlpha',0.3,'EdgeColor','k','LineWidth',0.5);
for k=1:size(FigureData.ExtendedEvents.BarPlots.CBVAvg,1)
    barHandle(k).CData=FigureData.ColorMaps.cmap((k+1),:);
end

figure(32); hold on
barHandle=bar(FigureData.ExtendedEvents.BarPlots.GCaMPAvg','FaceColor','flat','FaceAlpha',0.3,'EdgeColor','k','LineWidth',0.5);
for k=1:size(FigureData.ExtendedEvents.BarPlots.CBVAvg,1)
    barHandle(k).CData=FigureData.ColorMaps.cmap((k+1),:);
end

lgndTxt={'10s locomotion','15s locomotion','30s locomotion'};
for dataNum=1:size(dataTypes,1)
    
    figure(31);
    xVals=[];
    for k=1:size(FigureData.ExtendedEvents.BarPlots.(dataTypes{dataNum}).CBVVals,1)
        xVals(1:size(FigureData.ExtendedEvents.BarPlots.(dataTypes{dataNum}).CBVVals,2))=k;
        scatter((xVals+xOffset(dataNum)),FigureData.ExtendedEvents.BarPlots.(dataTypes{dataNum}).CBVVals(k,:),MarkerSize,'o','filled','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',FigureData.ColorMaps.cmap((dataNum+1),:));
    end
    
    if dataNum==size(dataTypes,1)
        xlim([0 (size(dataTypes,1)+1)]); xticks(1:1:size(dataTypes,1)); xticklabels(lgndTxt);
        ylabel('Z-units');
        title('Locomotion evoked blood volume change ');
%         savefig(['Population_Long_blood_volume_' theDate]);
%         saveas(gcf,['Population_Long_blood_volume_' theDate],'epsc');
    end
    
    
    figure(32);
    xVals=[];
    for k=1:size(FigureData.ExtendedEvents.BarPlots.(dataTypes{dataNum}).CBVVals,1)
        xVals(1:size(FigureData.ExtendedEvents.BarPlots.(dataTypes{dataNum}).CBVVals,2))=k;
        scatter((xVals+xOffset(dataNum)),FigureData.ExtendedEvents.BarPlots.(dataTypes{dataNum}).GCaMPVals(k,:),MarkerSize,'o','filled','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',FigureData.ColorMaps.cmap((dataNum+1),:));
    end
    
    if dataNum==size(dataTypes,1)
        xlim([0 (size(dataTypes,1)+1)]); xticks(1:1:size(dataTypes,1)); xticklabels(lgndTxt);
        ylabel('Z-units');
        title('Locomotion evoked GCaMP6s change ');
%         savefig(['Population_Long_GCaMP_' theDate]);
%         saveas(gcf,['Population_Long_GCaMP_' theDate],'epsc');
    end
    
    
    %% Fifteen second events same plot

    if dataNum==1
        lgndVals=[];
    end
    figure(33); hold on;
    colororder({'g','m'});
    yyaxis left;plot(FigureData.ExtendedEvents.FifteenSecEventAvg.(dataTypes{dataNum}).PlotTime,FigureData.ExtendedEvents.FifteenSecEventAvg.(dataTypes{dataNum}).GCaMP,'Color',thegreens{dataNum},'LineStyle','-','Marker','none');
    ylabel('Z-units');
    yyaxis right;plot(FigureData.ExtendedEvents.FifteenSecEventAvg.(dataTypes{dataNum}).PlotTime,FigureData.ExtendedEvents.FifteenSecEventAvg.(dataTypes{dataNum}).CBV,'Color',thereds{dataNum},'LineStyle','-','Marker','none');
    xlim([-5 15]); xlabel('Time(sec)'); ylabel('Z-units');title(' Population Average fifteen second locomotion evoked change ');
    TempTxt=strrep(dataTypes{dataNum},'_',' ');
    CBVTxt{1,(dataNum)}=[TempTxt ' blood volume'];
    GCaMPTxt{1,(dataNum)}=[TempTxt ' GCaMP'];
    
    if dataNum==length(dataTypes)
        lgndVals=[lgndVals,GCaMPTxt,CBVTxt];
        legend(lgndVals);
%         savefig(gcf,['Population_FifteenSecond_CBV_GFP_Evoked_Plots_' theDate]);
%         saveas(gcf,['Population_FifteenSecond_CBV_GFP_Evoked_Plots_' theDate],'epsc');
    end
    
    %% Ten second events same plot
    if dataNum==1
        lgndVals=[];
    end
    figure(34); hold on;
    colororder({'g','m'});
    yyaxis left;plot(FigureData.ExtendedEvents.TenSecEventAvg.(dataTypes{dataNum}).PlotTime,FigureData.ExtendedEvents.TenSecEventAvg.(dataTypes{dataNum}).GCaMP,'Color',thegreens{dataNum},'LineStyle','-','Marker','none');
    ylabel('Z-units');
    yyaxis right;plot(FigureData.ExtendedEvents.TenSecEventAvg.(dataTypes{dataNum}).PlotTime,FigureData.ExtendedEvents.TenSecEventAvg.(dataTypes{dataNum}).CBV,'Color',thereds{dataNum},'LineStyle','-','Marker','none');
    xlim([-5 10]); xlabel('Time(sec)'); ylabel('Z-units');title(' Population Average ten second locomotion evoked change ');
    TempTxt=strrep(dataTypes{dataNum},'_',' ');
    CBVTxt{1,(dataNum)}=[TempTxt ' blood volume'];
    GCaMPTxt{1,(dataNum)}=[TempTxt ' GCaMP'];
    
    if dataNum==length(dataTypes)
        lgndVals=[lgndVals,GCaMPTxt,CBVTxt];
        legend(lgndVals);
%         savefig(gcf,['Population_TenSecond_CBV_GFP_Evoked_Plots_' theDate]);
%         saveas(gcf,['Population_TenSecond_CBV_GFP_Evoked_Plots_' theDate],'epsc');
    end
    
    
    %% Thirty second events same plot
%     if dataNum==1
%         lgndVals=[];
%     end
%     figure(35); hold on;
%     colororder({'g','m'});
%     yyaxis left;plot(FigureData.ExtendedEvents.ThirtySecEventAvg.(dataTypes{dataNum}).PlotTime,FigureData.ExtendedEvents.ThirtySecEventAvg.(dataTypes{dataNum}).GCaMP,'Color',thegreens{dataNum},'LineStyle','-','Marker','none');
%     ylabel('Z-units');
%     yyaxis right;plot(FigureData.ExtendedEvents.ThirtySecEventAvg.(dataTypes{dataNum}).PlotTime,FigureData.ExtendedEvents.ThirtySecEventAvg.(dataTypes{dataNum}).CBV,'Color',thereds{dataNum},'LineStyle','-','Marker','none');
%     xlim([-5 30]); xlabel('Time(sec)'); ylabel('Z-units');title(' Population Average thirty second locomotion evoked change ');
%     TempTxt=strrep(dataTypes{dataNum},'_',' ');
%     CBVTxt{1,(dataNum)}=[TempTxt ' blood volume'];
%     GCaMPTxt{1,(dataNum)}=[TempTxt ' GCaMP'];
%     
%     if dataNum==length(dataTypes)
%         lgndVals=[lgndVals,GCaMPTxt,CBVTxt];
%         legend(lgndVals);
% %         savefig(gcf,['Population_ThirtySecond_CBV_GFP_Evoked_Plots_' theDate]);
% %         saveas(gcf,['Population_ThirtySecond_CBV_GFP_Evoked_Plots_' theDate],'epsc');
%     end
    clear xVals
end
end
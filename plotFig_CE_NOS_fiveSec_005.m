function plotFig_CE_NOS_fiveSec_005(FigureData)
dataTypes=fieldnames(FigureData.FiveSecFig.ResponseVol);
MarkerSize=100;
theDate=datetime;
theDate=char(theDate);
theDate=strrep(theDate,'-','');
theDate=strrep(theDate,':','');
theDate=strrep(theDate,' ','_');

for dataNum=1:size(dataTypes,1)
    titleTxt=strrep(dataTypes{dataNum},'_',' ');
    lgndTxt{dataNum}=strrep(dataTypes{dataNum},'_',' ');
    %% Response Volume    
    xVals(1:size(FigureData.FiveSecFig.ResponseVol.(dataTypes{dataNum}).SingleAnimalAvg.GCaMP,2))=dataNum; 
    figure(1); hold on
    bar(xVals(1),FigureData.FiveSecFig.ResponseVol.(dataTypes{dataNum}).PopulationAvg.CBV,'FaceColor',FigureData.ColorMaps.cmap((dataNum+1),:),'FaceAlpha',0.3,'EdgeColor','k','LineWidth',1);
    scatter(xVals,FigureData.FiveSecFig.ResponseVol.(dataTypes{dataNum}).SingleAnimalAvg.CBV,MarkerSize,'o','filled','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',FigureData.ColorMaps.cmap((dataNum+1),:));
    xlim([0 (size(dataTypes,1)+1)]); xticks(1:1:size(dataTypes,1)); xticklabels(lgndTxt);
    ylim([0 1.5]);
    ylabel('Z-units');
    title('Five second locomotion evoked blood volume change '); 
    if dataNum==size(dataTypes,1)
%     savefig(gcf,['Population_blood_volume_' theDate]);
%     saveas(gcf,['Population_blood_volume_' theDate],'epsc');
    end

    figure(2); hold on
    bar(xVals(1),FigureData.FiveSecFig.ResponseVol.(dataTypes{dataNum}).PopulationAvg.GCaMP,'FaceColor',FigureData.ColorMaps.cmap((dataNum+1),:),'FaceAlpha',0.3,'EdgeColor','k','LineWidth',1);
    scatter(xVals,FigureData.FiveSecFig.ResponseVol.(dataTypes{dataNum}).SingleAnimalAvg.GCaMP,MarkerSize,'o','filled','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',FigureData.ColorMaps.cmap((dataNum+1),:));
    xlim([0 (size(dataTypes,1)+1)]); xticks(1:1:size(dataTypes,1)); xticklabels(lgndTxt);
    ylim([0 3.5]);
    ylabel('Z-units');
    title('Five second locomotion evoked GCaMP6s change '); 
    if dataNum==size(dataTypes,1)
%     savefig(gcf,['Population_GCaMP_' theDate]);
%     saveas(gcf,['Population_GCaMP_' theDate],'epsc');
    end

    %% Five second events same plot
    if dataNum==1
        lgndVals=[];
    end
    thegreens={'#1A9850'; '#66BD63'; '#006837'};
    thereds={'#B30000'; '#D7301F'; '#7F0000'};
    figure(3); hold on;
    colororder({'g','m'});
    yyaxis left;plot(FigureData.FiveSecFig.LocomotionResponse.(dataTypes{dataNum}).PlotTime,...
        FigureData.FiveSecFig.LocomotionResponse.(dataTypes{dataNum}).GCaMP,'Color',thegreens{dataNum},'LineStyle','-','Marker','none');%FigureData.ColorMaps.cmap((dataNum+1),:)
    ylabel('Z-units');
    
    yyaxis right;plot(FigureData.FiveSecFig.LocomotionResponse.(dataTypes{dataNum}).PlotTime,...
        FigureData.FiveSecFig.LocomotionResponse.(dataTypes{dataNum}).CBV,'Color',thereds{dataNum},'LineStyle','-','Marker','none');%FigureData.ColorMaps.cmap((dataNum+1),:)
    ylabel('Z-units');
    
    xlim([-2 5]); xlabel('Time(sec)'); ylabel('Z-score');title(' Population Average five second locomotion evoked change ');
    TempTxt=strrep(dataTypes{dataNum},'_',' ');
    CBVTxt{1,dataNum}=[TempTxt ' blood volume'];
    GCaMPTxt{1,dataNum}=[TempTxt ' GCaMP'];
    if dataNum==length(dataTypes)
        lgndVals=[lgndVals,GCaMPTxt,CBVTxt];
        legend(lgndVals);
%         savefig(gcf,['Population_FiveSecond_CBV_GFP_Evoked_Plots_' theDate]);
%         saveas(gcf,['Population_FiveSecond_CBV_GFP_Evoked_Plots_' theDate],'epsc');
    end 
    clear xVals
end

%% Single Animal hSyn 1500um
figure(4);
ax1=subplot(5,1,1);
plot(FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.Trial.PlotTime,-100*FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.Trial.BallAcc,'k');%(1:end-1)
ylabel('Ball Velocity (cm/s)');
ylim([-1 8])
colororder({'g','m'});
ax2=subplot(5,1,(2:5)); hold on;
yyaxis right; plot(FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.Trial.PlotTime,FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.Trial.CBV,'m');ylim([-2 2]);ylabel('Z-units');
yyaxis left; plot(FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.Trial.PlotTime,FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.Trial.GCaMP,'g');ylim([-3 4]);
legend({'hSyn-GCaMP6s','Blood volume',},'Location','southwest');
ylabel('Z-units');

xticks(0:60:length(FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.Trial.GCaMP));
xlabel('Time (sec)');
linkaxes([ax1,ax2],'x');
xlim([1980 2100]);
title('5min period of behavior hSyn-GCaMP6s 1500um');
% savefig(gcf,['Single_Animal_hSyn_Trial_' theDate]);
% saveas(gcf,['Single_Animal_hSyn_Trial_' theDate],'epsc');


figure(5);hold on; 
colororder({'g','m'});
yyaxis left; plot(FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.PlotTime,FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.GCaMP,'Color','g','LineStyle','-','Marker','none'); ylim([-1 4]);
plot(FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.PlotTime,(FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.GCaMP+FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.GCaMPStd),'Color','g','LineStyle','--','Marker','none');
plot(FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.PlotTime,(FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.GCaMP-FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.GCaMPStd),'Color','g','LineStyle','--','Marker','none');
ylabel('Z-units');

yyaxis right;plot(FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.PlotTime,FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.CBV,'Color','m','LineStyle','-','Marker','none'); ylim([-0.5 1.5]);
plot(FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.PlotTime,(FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.CBV+FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.CBVStd),'Color','m','LineStyle','--','Marker','none');
plot(FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.PlotTime,(FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.CBV-FigureData.FiveSecFig.SingleAnimalPlots.Hsyn.AnAvg.CBVStd),'Color','m','LineStyle','--','Marker','none');
ylabel('Z-units');

xlabel('Time (sec)'); ylabel('Z-score'); title('Five second locomotion evoked change');xlim([-2 5]);
legend({'hSyn-GCaMP6s','GCaMP+std','GCaMP-std','Blood volume','Blood volume +std','Blood volume-std'});
% savefig(gcf,['Single_Animal_hSyn_fivesecAvg_' theDate]);
% saveas(gcf,['Single_Animal_hSyn_fivesecAvg_' theDate],'epsc');

%% Single Animal CaMKII 1500um

figure(6);
ax1=subplot(5,1,1);
plot(FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.Trial.PlotTime,-100*FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.Trial.BallAcc,'k');%(1:end-1)
ylabel('Ball Velocity (cm/s)');
ylim([-1 8])
ax2=subplot(5,1,(2:5)); hold on;
colororder({'g','m'});
yyaxis right; plot(FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.Trial.PlotTime,FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.Trial.CBV,'m');ylim([-2 2]);ylabel('Z-units');
ylim([-1 2]);
yyaxis left; plot(FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.Trial.PlotTime,FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.Trial.GCaMP,'g');ylim([-3 4]);ylabel('Z-units');
legend({'CaMKII-GCaMP6s','Blood volume'});
xticks(0:60:length(FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.Trial.GCaMP));
ylim([-3 3]);
xlabel('Time (sec)');
linkaxes([ax1,ax2],'x');
xlim([2520 2640]);
title('5min period of behavior CaMKII-GCaMP6s 1500um');
% savefig(gcf,['Single_Animal_CaMKII_Trial_' theDate]);
% saveas(gcf,['Single_Animal_CaMKII_Trial_' theDate],'epsc');

figure(7);hold on;
colororder({'g','m'});
yyaxis left; plot(FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.PlotTime,FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.GCaMP,'Color','g','LineStyle','-','Marker','none');
plot(FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.PlotTime,(FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.GCaMP+FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.GCaMPStd),'Color','g','LineStyle','--','Marker','none');
plot(FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.PlotTime,(FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.GCaMP-FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.GCaMPStd),'Color','g','LineStyle','--','Marker','none');
ylabel('Z-units');

yyaxis right; plot(FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.PlotTime,FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.CBV,'Color','m','LineStyle','-','Marker','none');
plot(FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.PlotTime,(FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.CBV+FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.CBVStd),'Color','m','LineStyle','--','Marker','none');
plot(FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.PlotTime,(FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.CBV-FigureData.FiveSecFig.SingleAnimalPlots.CaMKII.AnAvg.CBVStd),'Color','m','LineStyle','--','Marker','none');
ylabel('Z-units');

xlabel('Time (sec)'); ylabel('Z-score'); title('Five second locomotion evoked change');xlim([-2 5]);
legend({'CaMKII-GCaMP','GCaMP+std','GCaMP-std','Blood volume','Blood volume +std','Blood volume-std'},'Location','southeast');
% savefig(gcf,['Single_Animal_CaMKII_fivesecAvg_' theDate]);
% saveas(gcf,['Single_Animal_CaMKII_fivesecAvg_' theDate],'epsc');

%% Single Animal nNOS 1500um
figure(8);
ax1=subplot(5,1,1);
plot(FigureData.FiveSecFig.SingleAnimalPlots.nNOS.Trial.PlotTime,-100*FigureData.FiveSecFig.SingleAnimalPlots.nNOS.Trial.BallAcc,'k');%(1:end-1)
ylabel('Ball Velocity (cm/s)');
ylim([-10 20])
ax2=subplot(5,1,(2:5)); hold on;
colororder({'g','m'});
yyaxis right; plot(FigureData.FiveSecFig.SingleAnimalPlots.nNOS.Trial.PlotTime,FigureData.FiveSecFig.SingleAnimalPlots.nNOS.Trial.CBV,'m');ylim([-2 2]); ylabel('Z-units');
yyaxis left;plot(FigureData.FiveSecFig.SingleAnimalPlots.nNOS.Trial.PlotTime,FigureData.FiveSecFig.SingleAnimalPlots.nNOS.Trial.GCaMP,'g');ylim([-4 4]);ylabel('Z-units');
legend({'nNOS-GCaMP6s','Blood volume'},'Location','southeast');
ylabel('Z units');
xticks(0:60:length(FigureData.FiveSecFig.SingleAnimalPlots.nNOS.Trial.GCaMP));
xlabel('Time (sec)');
linkaxes([ax1,ax2],'x');
xlim([5940 6060]);
title('5min period of behavior nNos-GCaMP6s 1500um');
% savefig(gcf,['Single_Animal_nNos1500_Trial_' theDate]);
% saveas(gcf,['Single_Animal_nNos1500_Trial_' theDate],'epsc');


figure(9);hold on;
colororder({'g','m'});
yyaxis left; plot(FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.PlotTime,FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.GCaMP,'Color','g','LineStyle','-','Marker','none');
plot(FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.PlotTime,(FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.GCaMP+FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.GCaMPStd),'Color','g','LineStyle','--','Marker','none');
plot(FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.PlotTime,(FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.GCaMP-FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.GCaMPStd),'Color','g','LineStyle','--','Marker','none');
ylabel('Z-units');

yyaxis right; plot(FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.PlotTime,FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.CBV,'Color','m','LineStyle','-','Marker','none'); ylim([-0.5 1.5]);
plot(FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.PlotTime,(FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.CBV+FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.CBVStd),'Color','m','LineStyle','--','Marker','none');
plot(FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.PlotTime,(FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.CBV-FigureData.FiveSecFig.SingleAnimalPlots.nNOS.AnAvg.CBVStd),'Color','m','LineStyle','--','Marker','none');
ylabel('Z-units');

xlabel('Time (sec)'); ylabel('Z-score'); title('Five second locomotion evoked change'); xlim([-2 5]);
legend({'nNos-GCaMP','GCaMP+std','GCaMP-std','Blood volume','Blood volume +std','Blood volume-std'});
% savefig(gcf,['Single_Animal_nNos1500_fivesecAvg_' theDate]);
% saveas(gcf,['Single_Animal_nNos1500_fivesecAvg_' theDate],'epsc');
end

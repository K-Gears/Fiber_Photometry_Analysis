function plotFig_CE_NOS_SignalCorrection_005(FigureData)


%% PhotoBleaching correction
figure(41);plot(FigureData.SignalCorrection.CorrectMethod.PlotTime,FigureData.SignalCorrection.CorrectMethod.RawCBV,'m');hold on;  plot(FigureData.SignalCorrection.CorrectMethod.PlotTime,FigureData.SignalCorrection.CorrectMethod.PredictedCBVDecay,'c');title('Exponential fit of TRITC brightness'); xlabel('Time (sec)'); 
legend({'Raw TRITC brightness','Exponential Fit TRITC'}); xlim([0 FigureData.SignalCorrection.CorrectMethod.PlotTime(end)]);
% savefig('CAG_GFP_RAW_TRITC_FIT');
% saveas(gcf,'CAG_GFP_RAW_TRITC_FIT','epsc');

figure(42);plot(FigureData.SignalCorrection.CorrectMethod.PlotTime,FigureData.SignalCorrection.CorrectMethod.CorrectedCBV,'m'); title('Correction of TRITC brightness'); xlabel('Time (sec)'); 
legend('Corrected TRITC brightness'); xlim([0 FigureData.SignalCorrection.CorrectMethod.PlotTime(end)]);
% savefig('CAG_GFP_Corrected_TRITC');
%     saveas(gcf,'CAG_GFP_Corrected_TRITC','epsc');

figure(43);plot(FigureData.SignalCorrection.CorrectMethod.PlotTime,FigureData.SignalCorrection.CorrectMethod.RawGCaMP,'g');hold on;  plot(FigureData.SignalCorrection.CorrectMethod.PlotTime,FigureData.SignalCorrection.CorrectMethod.PredictedGCaMPDecay,'c');title('Exponential fit of GCaMP brightness'); xlabel('Time (sec)'); 
legend({'Raw GCaMP brightness','Exponential Fit GCaMP'}); xlim([0 FigureData.SignalCorrection.CorrectMethod.PlotTime(end)]);
% savefig('CAG_GFP_RAW_GCaMP_FIT');
% saveas(gcf,'CAG_GFP_RAW_GCaMP_FIT','epsc');

figure(44);plot(FigureData.SignalCorrection.CorrectMethod.PlotTime,FigureData.SignalCorrection.CorrectMethod.CorrectedGCaMP,'g');title('Correction of GCaMP brightness'); xlabel('Time (sec)'); 
legend('Corrected GCaMP brightness'); xlim([0 FigureData.SignalCorrection.CorrectMethod.PlotTime(end)]);
% savefig('CAG_GFP_Corrected_GCaMP');
% saveas(gcf,'CAG_GFP_Corrected_GCaMP','epsc');

LowPassData=[FigureData.SignalCorrection.CorrectMethod.CorrectedGCaMP,FigureData.SignalCorrection.CorrectMethod.CorrectedCBV];

figure(45); hold on;
yyaxis left; plot(FigureData.SignalCorrection.CorrectMethod.PlotTime,FigureData.SignalCorrection.CorrectMethod.CorrectedGCaMP,'g');
yyaxis right; plot(FigureData.SignalCorrection.CorrectMethod.PlotTime,FigureData.SignalCorrection.CorrectMethod.CorrectedCBV,'m');
title('Detrend Fiber Data'); xlabel('Time (sec)'); ylabel('Fluorescence (A.U.)'); 
legend({'eGFP','TRITC'}); xlim([0 FigureData.SignalCorrection.CorrectMethod.PlotTime(end)]);
% savefig('CAG_GFP_Rescaled_Uncorrected_Data');

%% Correction of Hemodynamic GFP attenuation
figure(46); hold on;
imagesc(FigureData.SignalCorrection.CorrectionHist.BinEdges,FigureData.SignalCorrection.CorrectionHist.BinEdges,FigureData.SignalCorrection.CorrectionHist.BinCounts');
caxis([0 10]);
h=colorbar('eastoutside');
axis xy;
ylabel('Observed GFP');
xlabel('Observed TRITC');
title('Normalized histogram of CBV vs GFP');
xlim([0 1]);
ylim([0 1]);
plot(FigureData.SignalCorrection.CorrectionHist.BinEdges,FigureData.SignalCorrection.CorrectionHist.LinearFit,'c','LineWidth',2);
legend({'Linear fit, peak column values'});
set(get(h,'label'),'string','Percentage of all column counts in bin (%)');
%% Single Animal eGFP 1500um

figure(47);
ax1=subplot(5,1,1);
plot(FigureData.SignalCorrection.eGFPTrial.PlotTime,100*abs(FigureData.SignalCorrection.eGFPTrial.BallAcc),'k');%(1:end-1)
ylim([0 1.5])
ax2=subplot(5,1,(2:5)); hold on;
colororder({'g','m'});
yyaxis right; plot(FigureData.SignalCorrection.eGFPTrial.PlotTime,FigureData.SignalCorrection.eGFPTrial.CorrectedCBV,'m'); ylim([-1.5 3]); ylabel('Z-units');
yyaxis left; plot(FigureData.SignalCorrection.eGFPTrial.PlotTime,FigureData.SignalCorrection.eGFPTrial.CorrectedGCaMP,'g');
legend({'TRITC','eGFP'});
ylabel('Z-units');
xticks(0:10:length(FigureData.SignalCorrection.eGFPTrial.CorrectedGCaMP));
ylim([-2 2]);
xlabel('Time (sec)');
linkaxes([ax1,ax2],'x');
xlim([680 740]);
title('5min period of behavior CAG-eGFP corrected 1500um');
% savefig(gcf,['Single_Animal_eGFP_Trial_' theDate]);
% saveas(gcf,['Single_Animal_eGFP_Trial_' theDate],'epsc');

figure(48);
ax1=subplot(5,1,1);
plot(FigureData.SignalCorrection.eGFPTrial.PlotTime,100*abs(FigureData.SignalCorrection.eGFPTrial.BallAcc),'k');%(1:end-1)
ylim([0 1.5])
ax2=subplot(5,1,(2:5)); hold on;
colororder({'g','m'});
yyaxis right; plot(FigureData.SignalCorrection.eGFPTrial.PlotTime,FigureData.SignalCorrection.eGFPTrial.CorrectedCBV,'m'); ylim([-1.5 3]); ylabel('Z-units');
yyaxis left; plot(FigureData.SignalCorrection.eGFPTrial.PlotTime,FigureData.SignalCorrection.eGFPTrial.UncorrectedGCaMP,'g');
legend({'TRITC','eGFP'});
ylabel('Z-units');
xticks(0:10:length(FigureData.SignalCorrection.eGFPTrial.UncorrectedGCaMP));
ylim([-2 2]);
xlabel('Time (sec)');
linkaxes([ax1,ax2],'x');
xlim([680 740]);
title('5min period of behavior CAG-eGFP raw 1500um');
% savefig(gcf,['Single_Animal_eGFP_Trial_Raw' theDate]);
% saveas(gcf,['Single_Animal_eGFP_Trial_Raw' theDate],'epsc');




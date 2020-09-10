function fiberphotometry_figures_004(FigureData)

if ischar(FigureData)
    load(FigureData);
end
plotFig_CE_NOS_fiveSec_005(FigureData);
plotFig_CE_NOS_Xcorr_005(FigureData);
plotFig_CE_NOS_LongEvents_005(FigureData);
plotFig_CE_NOS_SignalCorrection_005(FigureData);

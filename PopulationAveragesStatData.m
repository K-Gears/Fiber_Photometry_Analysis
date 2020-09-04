function [PopStats]=PopulationAveragesStatData(FigureData)
%% GCaMP Evoked Intensity
PopStats.CaMKII.GCaMP5sAvg=mean(FigureData.Stats.StatsData.GCaMPPeak5s(1:5));
PopStats.CaMKII.GCaMP5sStd=std(FigureData.Stats.StatsData.GCaMPPeak5s(1:5));

PopStats.CaMKII.GCaMP10sAvg=mean(FigureData.Stats.StatsData.GCaMPPeak10s(1:5));
PopStats.CaMKII.GCaMP10sStd=std(FigureData.Stats.StatsData.GCaMPPeak10s(1:5));

PopStats.CaMKII.GCaMP15sAvg=mean(FigureData.Stats.StatsData.GCaMPPeak15s(1:5));
PopStats.CaMKII.GCaMP15sStd=std(FigureData.Stats.StatsData.GCaMPPeak15s(1:5));

PopStats.CaMKII.GCaMP30sAvg=mean(FigureData.Stats.StatsData.GCaMPPeak30s(1:5));
PopStats.CaMKII.GCaMP30sStd=std(FigureData.Stats.StatsData.GCaMPPeak30s(1:5));

PopStats.NOS.GCaMP5sAvg=mean(FigureData.Stats.StatsData.GCaMPPeak5s(6:14));
PopStats.NOS.GCaMP5sStd=std(FigureData.Stats.StatsData.GCaMPPeak5s(6:14));

PopStats.NOS.GCaMP10sAvg=mean(FigureData.Stats.StatsData.GCaMPPeak10s(6:14));
PopStats.NOS.GCaMP10sStd=std(FigureData.Stats.StatsData.GCaMPPeak10s(6:14));

PopStats.NOS.GCaMP15sAvg=mean(FigureData.Stats.StatsData.GCaMPPeak15s(6:14));
PopStats.NOS.GCaMP15sStd=std(FigureData.Stats.StatsData.GCaMPPeak15s(6:14));

PopStats.NOS.GCaMP30sAvg=mean(FigureData.Stats.StatsData.GCaMPPeak30s(6:14));
PopStats.NOS.GCaMP30sStd=std(FigureData.Stats.StatsData.GCaMPPeak30s(6:14));

PopStats.hSyn.GCaMP5sAvg=mean(FigureData.Stats.StatsData.GCaMPPeak5s(15:18));
PopStats.hSyn.GCaMP5sStd=std(FigureData.Stats.StatsData.GCaMPPeak5s(15:18));

PopStats.hSyn.GCaMP10sAvg=mean(FigureData.Stats.StatsData.GCaMPPeak10s(15:18));
PopStats.hSyn.GCaMP10sStd=std(FigureData.Stats.StatsData.GCaMPPeak10s(15:18));

PopStats.hSyn.GCaMP15sAvg=mean(FigureData.Stats.StatsData.GCaMPPeak15s(15:18));
PopStats.hSyn.GCaMP15sStd=std(FigureData.Stats.StatsData.GCaMPPeak15s(15:18));

PopStats.hSyn.GCaMP30sAvg=mean(FigureData.Stats.StatsData.GCaMPPeak30s(15:18));
PopStats.hSyn.GCaMP30sStd=std(FigureData.Stats.StatsData.GCaMPPeak30s(15:18));

%% TRITC Evoked Intensity
PopStats.CaMKII.CBV5sAvg=mean(FigureData.Stats.StatsData.CBVPeak5s(1:5));
PopStats.CaMKII.CBV5sStd=std(FigureData.Stats.StatsData.CBVPeak5s(1:5));

PopStats.CaMKII.CBV10sAvg=mean(FigureData.Stats.StatsData.CBVPeak10s(1:5));
PopStats.CaMKII.CBV10sStd=std(FigureData.Stats.StatsData.CBVPeak10s(1:5));

PopStats.CaMKII.CBV15sAvg=mean(FigureData.Stats.StatsData.CBVPeak15s(1:5));
PopStats.CaMKII.CBV15sStd=std(FigureData.Stats.StatsData.CBVPeak15s(1:5));

PopStats.CaMKII.CBV30sAvg=mean(FigureData.Stats.StatsData.CBVPeak30s(1:5));
PopStats.CaMKII.CBV30sStd=std(FigureData.Stats.StatsData.CBVPeak30s(1:5));

PopStats.NOS.CBV5sAvg=mean(FigureData.Stats.StatsData.CBVPeak5s(6:14));
PopStats.NOS.CBV5sStd=std(FigureData.Stats.StatsData.CBVPeak5s(6:14));

PopStats.NOS.CBV10sAvg=mean(FigureData.Stats.StatsData.CBVPeak10s(6:14));
PopStats.NOS.CBV10sStd=std(FigureData.Stats.StatsData.CBVPeak10s(6:14));

PopStats.NOS.CBV15sAvg=mean(FigureData.Stats.StatsData.CBVPeak15s(6:14));
PopStats.NOS.CBV15sStd=std(FigureData.Stats.StatsData.CBVPeak15s(6:14));

PopStats.NOS.CBV30sAvg=mean(FigureData.Stats.StatsData.CBVPeak30s(6:14));
PopStats.NOS.CBV30sStd=std(FigureData.Stats.StatsData.CBVPeak30s(6:14));

PopStats.hSyn.CBV5sAvg=mean(FigureData.Stats.StatsData.CBVPeak5s(15:18));
PopStats.hSyn.CBV5sStd=std(FigureData.Stats.StatsData.CBVPeak5s(15:18));

PopStats.hSyn.CBV10sAvg=mean(FigureData.Stats.StatsData.CBVPeak10s(15:18));
PopStats.hSyn.CBV10sStd=std(FigureData.Stats.StatsData.CBVPeak10s(15:18));

PopStats.hSyn.CBV15sAvg=mean(FigureData.Stats.StatsData.CBVPeak15s(15:18));
PopStats.hSyn.CBV15sStd=std(FigureData.Stats.StatsData.CBVPeak15s(15:18));

PopStats.hSyn.CBV30sAvg=mean(FigureData.Stats.StatsData.CBVPeak30s(15:18));
PopStats.hSyn.CBV30sStd=std(FigureData.Stats.StatsData.CBVPeak30s(15:18));

%% Xcorr Peaks/Lags
PopStats.CaMKII.CorrCoeffAvg=mean(FigureData.Stats.StatsData.XCorrCoeffPeak(1:5));
PopStats.CaMKII.CorrCoeffStd=std(FigureData.Stats.StatsData.XCorrCoeffPeak(1:5));

PopStats.CaMKII.CorrLagAvg=mean(FigureData.Stats.StatsData.XCorrCoeffLag(1:5));
PopStats.CaMKII.CorrLagStd=std(FigureData.Stats.StatsData.XCorrCoeffLag(1:5));

PopStats.NOS.CorrCoeffAvg=mean(FigureData.Stats.StatsData.XCorrCoeffPeak(6:14));
PopStats.NOS.CorrCoeffStd=std(FigureData.Stats.StatsData.XCorrCoeffPeak(6:14));

PopStats.NOS.CorrLagAvg=mean(FigureData.Stats.StatsData.XCorrCoeffLag(6:14));
PopStats.NOS.CorrLagStd=std(FigureData.Stats.StatsData.XCorrCoeffLag(6:14));

PopStats.hSyn.CorrCoeffAvg=mean(FigureData.Stats.StatsData.XCorrCoeffPeak(15:18));
PopStats.hSyn.CorrCoeffStd=std(FigureData.Stats.StatsData.XCorrCoeffPeak(15:18));

PopStats.hSyn.CorrLagAvg=mean(FigureData.Stats.StatsData.XCorrCoeffLag(15:18));
PopStats.hSyn.CorrLagStd=std(FigureData.Stats.StatsData.XCorrCoeffLag(15:18));

%% Raw Signal Intensity
PopStats.CaMKII.RawGCaMPAvg=mean(FigureData.Stats.StatsData.RawGCaMPAvg(1:5));
PopStats.CaMKII.RawGCaMPStd=std(FigureData.Stats.StatsData.RawGCaMPAvg(1:5));

PopStats.CaMKII.RawCBVAvg=mean(FigureData.Stats.StatsData.RawCBVAvg(1:5));
PopStats.CaMKII.RawCBVStd=std(FigureData.Stats.StatsData.RawCBVAvg(1:5));

PopStats.NOS.RawGCaMPAvg=mean(FigureData.Stats.StatsData.RawGCaMPAvg(6:14));
PopStats.NOS.RawGCaMPStd=std(FigureData.Stats.StatsData.RawGCaMPAvg(6:14));

PopStats.NOS.RawCBVAvg=mean(FigureData.Stats.StatsData.RawCBVAvg(6:14));
PopStats.NOS.RawCBVStd=std(FigureData.Stats.StatsData.RawCBVAvg(6:14));

PopStats.hSyn.RawGCaMPAvg=mean(FigureData.Stats.StatsData.RawGCaMPAvg(15:18));
PopStats.hSyn.RawGCaMPStd=std(FigureData.Stats.StatsData.RawGCaMPAvg(15:18));

PopStats.hSyn.RawCBVAvg=mean(FigureData.Stats.StatsData.RawCBVAvg(15:18));
PopStats.hSyn.RawCBVStd=std(FigureData.Stats.StatsData.RawCBVAvg(15:18));

%% Raw Signal Intensity STD
PopStats.CaMKII.RawStdGCaMPAvg=mean(FigureData.Stats.StatsData.RawGCaMPStd(1:5));
PopStats.CaMKII.RawStdGCaMPStd=std(FigureData.Stats.StatsData.RawGCaMPStd(1:5));

PopStats.CaMKII.RawStdCBVAvg=mean(FigureData.Stats.StatsData.RawCBVStd(1:5));
PopStats.CaMKII.RawStdCBVStd=std(FigureData.Stats.StatsData.RawCBVStd(1:5));

PopStats.NOS.RawStdGCaMPAvg=mean(FigureData.Stats.StatsData.RawGCaMPStd(6:14));
PopStats.NOS.RawStdGCaMPStd=std(FigureData.Stats.StatsData.RawGCaMPStd(6:14));

PopStats.NOS.RawStdCBVAvg=mean(FigureData.Stats.StatsData.RawCBVStd(6:14));
PopStats.NOS.RawStdCBVStd=std(FigureData.Stats.StatsData.RawCBVStd(6:14));

PopStats.hSyn.RawStdGCaMPAvg=mean(FigureData.Stats.StatsData.RawGCaMPStd(15:18));
PopStats.hSyn.RawStdGCaMPStd=std(FigureData.Stats.StatsData.RawGCaMPStd(15:18));

PopStats.hSyn.RawStdCBVAvg=mean(FigureData.Stats.StatsData.RawCBVStd(15:18));
PopStats.hSyn.RawStdCBVStd=std(FigureData.Stats.StatsData.RawCBVStd(15:18));
function [GLME_Output]=GLME_FiberData_Fixed(Locomotion_Data_Table)
% gcamp_formula='GCaMPPeak ~ 1+AnimalGenotype(1|AnimalProbeDepth)';
% cbv_formula='CBVPeak ~ 1+AnimalGenotype+(1|AnimalProbeDepth)';
% gcamp_stats=fitglme(Locomotion_Data_Table,gcamp_formula);%relative to Rest "Intercept"
% cbv_stats=fitglme(Locomotion_Data_Table,cbv_formula);%relative to Rest "Intercept"

%% Locomotion Evoked Statistics GCaMP
FiveSec_GCaMP6s_Stats=fitglme(Locomotion_Data_Table,'GCaMPPeak5s ~ 1+AnimalGenotype+AnimalProbeDepth')
TenSec_GCaMP6s_Stats=fitglme(Locomotion_Data_Table,'GCaMPPeak10s ~ 1+AnimalGenotype+AnimalProbeDepth')
FifteenSec_GCaMP6s_Stats=fitglme(Locomotion_Data_Table,'GCaMPPeak15s ~ 1+AnimalGenotype+AnimalProbeDepth')
ThirtySec_GCaMP6s_Stats=fitglme(Locomotion_Data_Table,'GCaMPPeak30s ~ 1+AnimalGenotype+AnimalProbeDepth')

%% Locomotion Evoked Statistics CBV
FiveSec_CBV_Stats=fitglme(Locomotion_Data_Table,'CBVPeak5s ~ 1+AnimalGenotype+AnimalProbeDepth')
TenSec_CBV_Stats=fitglme(Locomotion_Data_Table,'CBVPeak10s ~ 1+AnimalGenotype+AnimalProbeDepth')
FifteenSec_CBV_Stats=fitglme(Locomotion_Data_Table,'CBVPeak15s ~ 1+AnimalGenotype+AnimalProbeDepth')
ThirtySec_CBV_Stats=fitglme(Locomotion_Data_Table,'CBVPeak30s ~ 1+AnimalGenotype+AnimalProbeDepth')

%% Aggregate Locomotion Stats
GLME_Output.Locomotion.FiveSec.GCaMP6s_Stats=FiveSec_GCaMP6s_Stats;
GLME_Output.Locomotion.TenSec.GCaMP6s_Stats=TenSec_GCaMP6s_Stats;
GLME_Output.Locomotion.FifteenSec.GCaMP6s_Stats=FifteenSec_GCaMP6s_Stats;
GLME_Output.Locomotion.ThirtySec.GCaMP6s_Stats=ThirtySec_GCaMP6s_Stats;

GLME_Output.Locomotion.FiveSec.CBV_Stats=FiveSec_CBV_Stats;
GLME_Output.Locomotion.TenSec.CBV_Stats=TenSec_CBV_Stats;
GLME_Output.Locomotion.FifteenSec.CBV_Stats=FifteenSec_CBV_Stats;
GLME_Output.Locomotion.ThirtySec.CBV_Stats=ThirtySec_CBV_Stats;

%% Raw Signal Statistics GCaMP
RawAvg_GCaMP6s_Stats=fitglme(Locomotion_Data_Table,'RawGCaMPAvg ~ 1+AnimalGenotype+AnimalProbeDepth')
RawStd_GCaMP6s_Stats=fitglme(Locomotion_Data_Table,'RawGCaMPStd ~ 1+AnimalGenotype+AnimalProbeDepth','Distribution','Normal','Link','identity')

%% Raw Signal Statistics CBV
RawAvg_CBV_Stats=fitglme(Locomotion_Data_Table,'RawCBVAvg ~ 1+AnimalGenotype+AnimalProbeDepth')
RawStd_CBV_Stats=fitglme(Locomotion_Data_Table,'RawCBVStd ~ 1+AnimalGenotype+AnimalProbeDepth','Distribution','Normal','Link','identity')


%% Aggregate Raw Signal Statistics
GLME_Output.RawSignalNoise.RawAvg.GCaMP6s_Stats=RawAvg_GCaMP6s_Stats;
GLME_Output.RawSignalNoise.RawStd.GCaMP6s_Stats=RawStd_GCaMP6s_Stats;

GLME_Output.RawSignalNoise.RawAvg.CBV_Stats=RawAvg_CBV_Stats;
GLME_Output.RawSignalNoise.RawStd.CBV_Stats=RawStd_CBV_Stats;

%% XCorr Statistics
Peaks_Stats=fitglme(Locomotion_Data_Table,'XCorrCoeffPeak ~ 1+AnimalGenotype+AnimalProbeDepth')
Lags_Stats=fitglme(Locomotion_Data_Table,'XCorrCoeffLag ~ 1+AnimalGenotype+AnimalProbeDepth')


%% Aggregate XCorr Statistics
GLME_Output.Xcorr.Peaks_Stats=Peaks_Stats;
GLME_Output.Xcorr.Lags_Stats=Lags_Stats;

TenSec_GCaMP6s_Stats=fitglme(Locomotion_Data_Table((6:14),:),'GCaMPPeak10s ~ 1+GCaMPPeak5s+AnimalProbeDepth')%(1|AnimalProbeDepth)')
% [psi,dispersion,stats] = covarianceParameters(TenSec_GCaMP6s_Stats);

FifteenSec_GCaMP6s_Stats=fitglme(Locomotion_Data_Table((6:14),:),'GCaMPPeak15s ~ 1+GCaMPPeak5s+AnimalProbeDepth')%(1|AnimalProbeDepth)')
% [psi,dispersion,stats] = covarianceParameters(FifteenSec_GCaMP6s_Stats);

ThirtySec_GCaMP6s_Stats=fitglme(Locomotion_Data_Table((6:14),:),'GCaMPPeak30s ~ 1+GCaMPPeak5s+AnimalProbeDepth')%(1|AnimalProbeDepth)')

GLME_Output.Locomotion.TenSec.GCaMP_Duration_Stats=TenSec_GCaMP6s_Stats;
GLME_Output.Locomotion.FifteenSec.GCaMP_Duration_Stats=FifteenSec_GCaMP6s_Stats;
GLME_Output.Locomotion.ThirtySec.GCaMP_Duration_Stats=ThirtySec_GCaMP6s_Stats;
end
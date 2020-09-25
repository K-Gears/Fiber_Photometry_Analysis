function [correctionConstant,initCorrCoeffs,adjCorrCoeffs]= MinimizeCorrCoeff(GCaMP,CBV,Fs)
%Inputs:
% GCaMP: this is the time series of the GCaMP signal that has been
% corrected for photo-bleaching
% CBV: This is the time series of the TRITC signal that has been corrected
% for metabolic signal decay
% Fs: signal acquisition rate
%
%Outputs:
%correctionConstant: scaling factor to be applied to TRITC signal prior to
%addition to GCaMP signal

[z,p,k]=butter(3,1/(0.5*Fs),'low');
[sos,g]=zp2sos(z,p,k); %1Hz lowpass filter to apply to signals before analysis
CBV=filtfilt(sos,g,CBV);
GCaMP=filtfilt(sos,g,GCaMP);

plotTime=(1:length(CBV))/1200;
figure(101);
ax1=subplot(2,1,1);plot(plotTime,detrend(CBV),'Color','m','LineWidth',1);
hold on; plot(plotTime,detrend(GCaMP),'Color','g','LineWidth',1); xlabel('Time (sec)'); ylabel('a.u'); title('Uncorrected GCaMP and TRITC');
 legend({'TRITC','GCaMP'});
[initCorrCoeffs]=xcorr(detrend(CBV),detrend(GCaMP),0,'coeff');

%% Set up minimization to estimate correction constant
initConst=0.2; %Best guess of scaling factor for hemodynamic correct to start minimization

options = optimset('Display','off',...      % display no output
    'FunValCheck','off',...  % check objective values
    'MaxFunEvals', 1000,...  % max number of function evaluations allowed
    'MaxIter', 1000,...      % max number of iteration allowed
    'TolFun',1e-8,...        % termination tolerance on the function value
    'TolX',1e-8,...          % termination tolerance on x
    'UseParallel','always'); % always use parallel computation

correctionConstant=fminsearch(@(y) ModelCoeffMin(GCaMP,CBV,y),initConst,options);

adjustedGCaMP=detrend(GCaMP)+correctionConstant*detrend(CBV);

figure(101);
ax2=subplot(2,1,2);plot(plotTime,detrend(CBV),'Color','m','LineWidth',1);
hold on; plot(plotTime,detrend(adjustedGCaMP),'Color','g','LineWidth',1); 
xlabel('Time (sec)'); ylabel('a.u'); title('Corrected GCaMP and TRITC');
legend({'TRITC','GCaMP'});
linkaxes([ax1,ax2],'x'); xlim([0 7000]);

[adjCorrCoeffs]=xcorr(detrend(CBV),detrend(adjustedGCaMP),0,'coeff');
end

%% Function to Minimize
function corrCoeff=ModelCoeffMin(GCaMP,CBV,y)
% ModelCoeffMin -- This function seeks to minimize the correlation
% coefficient at zero lag between GCaMP signals and CBV signals
% Inputs:
% GCaMP: this is the time series of the GCaMP signal that has been
% corrected for photo-bleaching
% CBV: This is the time series of the TRITC signal that has been corrected
% for metabolic signal decay
% StartVal: Decimal estimate (0-1) for scaling factor to minimize corrCoeff at
% lag==0
%
%
% Outputs:
% corrCoeff: squared correlation coefficient of xcorr(CBV,correctedGCaMP,0,'coeff')
%
%Written By: Kyle Gheres (kwg5014@gmail.com) 09-17-20

theScale=y;

correctedGCaMP=detrend(GCaMP)+(theScale*detrend(CBV)); % attempt to remove reflected copy of TRITC signal

corrCoeff=xcorr(detrend(CBV),detrend(correctedGCaMP),0,'coeff')^2; %calculate and square zero lagged correlation coefficeint between correctGCaMP and TRITC signal

end
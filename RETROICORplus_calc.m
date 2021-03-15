%--------------------------------------------------------------------------
%RETROICORplus_calc
%
%Uses heart rate and respiration data to create nuisance regressors
%for physiological pulsations. 
%
%Output:
%CPR: cardiac phase regressors (2 (sin+cos) * order regressors)
%RPR: respiratory phase regressors (2 (sin+cos) * order regressors)
%NR: other nuisance regressors:
%   HRF: heart rate frequency at different temporal delays (1 regressors per time shift)
%   RVT: Frequency times amplitude of respiration at different temporal delays (1 regressors per time shift) 

%EJH 2010-21
%--------------------------------------------------------------------------
function [CPR,RPR,NR]=RETROICORplus_calc(TTLlines,Peaklines,Pulsedat,Respdat,sR,RETROICORplus_defaults)


%---------------------------------------------------------------
%Get settings from defaults
%---------------------------------------------------------------
HRFwinlen = RETROICORplus_defaults.HRFwinlen;   %Window length in sec for HRF calculation
RVTwinlen = RETROICORplus_defaults.RVTwinlen;   %Window length in sec for RVT calculations
Respphasewinlen = RETROICORplus_defaults.Respphasewinlen; %Window length in sec for respiratory phase estimation
fOrder = RETROICORplus_defaults.fOrder;         %Fourier order for retroicor regressors
TS_RVT = RETROICORplus_defaults.TS_RVT;         %Time shift(s) for RVT
TS_HRF = RETROICORplus_defaults.TS_HRF;         %Time shift(s) for HR
ndat = numel(Pulsedat);                         %Get number of raw data points

%---------------------------------------------------------------
%Calculate cardiac phase for each TR
%---------------------------------------------------------------

%Initialize empty vectors to fill up later
HRphasedat=zeros(1,ndat);             %Cardiac phase data channel
IBIdat=zeros(1,ndat);                 %IBI data channel

%First take all pulse peaks and calculate the phase at each time point
%Loop over pulses to calculate phase and IBI data
for i=2:length(Peaklines)
    HRphasedat(Peaklines(i-1):Peaklines(i)-1) = ...
        linspace(0,2*pi,diff(Peaklines(i-1:i)));
    IBIdat(Peaklines(i-1):Peaklines(i)-1) = diff(Peaklines(i-1:i));
end

%Extrapolate start and end windows of IBIdat
IBIdat(1:Peaklines(1))= IBIdat(Peaklines(1));
IBIdat(Peaklines(end):end)=IBIdat(Peaklines(end)-1);

%Get cardiac phase for each TTL pulse
%20170228: avoid crash if there are TTL pulses beyond the end of the file
%which can happen if the recording is stopped immediately after scanning,
%and the filters lead to a shortening of the physiological recording
%old: HRPhase_TTL = HRphasedat(TTLlines);

%20210311: bug in these lines: it is now feeding the scan trigger samples
%(HRPhase_TTL) themselves into the CPR fourier expansion. that should of course be the
%cardiac phase at those scan triggers, so indexing into the HRphasedat.
%Old lines:
%HRPhase_TTL = zeros(1,numel(TTLlines));
%HRPhase_TTL(1:sum(TTLlines<=numel(HRphasedat))) = ...
%    TTLlines(TTLlines<=numel(HRphasedat));
%corrected:
HRPhase_TTL = zeros(1,numel(TTLlines));
HRPhase_TTL(1:sum(TTLlines<=numel(HRphasedat))) = ...
    HRphasedat(TTLlines(TTLlines<=numel(HRphasedat)));



%---------------------------------------------------------------
%Calculate the respiration phase for each TR
%---------------------------------------------------------------

warning off;

sRespdat = (Respdat - min(Respdat)) ./ (max(Respdat)-min(Respdat));   %Rescale the respiration to 0 to 1
[H, b] = hist(sRespdat, 100);               %Histogram data with 100 bins +frequencies in sRespdat
Respphasewinsam = Respphasewinlen*sR;       %Window length in samples for phase estimation
for i=1:numel(TTLlines)                     %Loop over TTLs
    cTTLloc = TTLlines(i);                  %Get the index in the raw data of the current TTL
    WinSt=max(1,cTTLloc-ceil(Respphasewinsam/2));
    WinEn=min(ndat,cTTLloc+ceil(Respphasewinsam/2));    
    WinResp = Respdat(WinSt:WinEn);           %WinRESP is the raw Resp channel within the window
    
    if numel(WinResp)<(Respphasewinsam/2)   %If less than half a window of data is available
        RespSigns_TTL(i) = 0;               %Set phase and sign to zero if no data
        RespPhase_TTL(i) = 0;
    else %if enough data

        %Determine if inhale or exhale by estimating the sign of the slope
        %in the middle of the respiration window
        X=[ones(numel(WinResp),1),(WinSt:WinEn)',((WinSt:WinEn).^2)']; %Setup DM with constant, line + square     %%%% WARNING is close to singular or badly scaled MvB
        betas = inv(X'*X)*(X')*WinResp';          %Fit DM onto data
        RespSigns_TTL(i) = sign(cTTLloc*2*betas(3) + betas(2));  %Get sign of slope in the middle of the window
        f = find(b <= sRespdat(cTTLloc));
        temp = sum(H(f));       %how many datapoints are below the current value
        RespPhase_TTL(i) = (pi * RespSigns_TTL(i) * temp / length(sRespdat))+pi;
        
    end %if enough data
end %loop over TTLs

warning on;

%---------------------------------------------------------------
%Calculation of HRF with windows surrounding TRs + time shift(s)
%---------------------------------------------------------------
HRFwinsam=HRFwinlen*sR;                           %Window Length in samples for heart rate
for cTS_HRF=1:numel(TS_HRF)                       %Loop over any time shifts requested in defaults

    tsTTLlines = TTLlines - TS_HRF(cTS_HRF)*sR;   %Apply the default time shift for heart rate
    HRF_tsTTL=[];                               %Clear variable
    for i=1:numel(tsTTLlines)                   %Loop over scanner pulse TTLs
        ctsTTLloc = tsTTLlines(i);              %Get the index in the raw data of the current TTL
        WinSt=max(1,ctsTTLloc-ceil(HRFwinsam/2));
        WinEn=min(ndat,ctsTTLloc+ceil(HRFwinsam/2));    
        WinIBI = IBIdat(WinSt:WinEn);           %WinIBI is the IBIchannel within the window
        WinIBI(find(WinIBI==0)) = [];           %remove zeros (although these should not be there)
        if mean(WinIBI)>0
            HRF_tsTTL(i)= 1/(mean(WinIBI)/sR); %HRF: mean heart rate freq in Hz
        else
            HRF_tsTTL(i)=NaN;           %Put in NaN if there's no data
        end
    end

    %Extrapolate beginning and end if NaNs
    if sum(~isnan(HRF_tsTTL))==0                %If there's no HR data at all
        HRF_tsTTL(:)=0;                         %Set to all zero
    end
    if isnan(HRF_tsTTL(1))                                          %If there's NaNs at the beginning
        stDatind = find(~isnan(HRF_tsTTL));                         %Get the first real datapoint
        HRF_tsTTL(1:stDatind(1)) = HRF_tsTTL(stDatind(1));          %Extrapolate to beginning
    end
    if isnan(HRF_tsTTL(end))                                        %If there's NaNs at the end
        endDatind = find(~isnan(HRF_tsTTL));                        %Get the last real datapoint
        HRF_tsTTL(endDatind(end):end) = HRF_tsTTL(endDatind(end)); %Extrapolate to end
    end
    
    %Collect the time shifted HRF vectors
    HRF_allTS(:,cTS_HRF) = HRF_tsTTL';

end %Loop over HR time shifts
    
%---------------------------------------------------------------    
%Calculation of RVT with windows surrounding TRs + time shift(s)
%---------------------------------------------------------------
RVTwinsam = round(RVTwinlen * sR);
for cTS_RVT=1:numel(TS_RVT)                       %Loop over any time shifts requested in defaults

    tsTTLlines = TTLlines - TS_RVT(cTS_RVT)*sR;   %Apply the default time shift for heart rate
    RVT_tsTTL=[];                                 %Clear variables
    WinEn=[];                                      %added MvB                     
    WinSt=[];                                      %added MvB
    for i=1:numel(tsTTLlines)                     %Loop over scanner pulse TTLs
        ctsTTLloc = tsTTLlines(i);                %Get the index in the raw data of the current TTL
        WinSt=max(1,ctsTTLloc-ceil(RVTwinsam/2));
        WinEn=min(ndat,ctsTTLloc+ceil(RVTwinsam/2));
        WinResp = Respdat(WinSt:WinEn);           %WinRESP is the raw Resp channel within the window
        
        %If there is sufficient data within this window (at least half a
        %window length) do calculations, otherwise put in NaN
        if numel(WinResp)<(RVTwinsam/2)
            RVT_tsTTL(i) = NaN;
        else %Run the RVT calculation
            %Detrend the window
            P = polyfit(1:numel(WinResp),WinResp,1);
            WinResp=WinResp-((1:numel(WinResp)).*P(1)+P(2));
            %Calculate amplitude/frequency and RVT
            fft0 = fft(WinResp);                                %fourier transform respiration window
            fft0 = fft0(1:floor(length(fft0) / 2));
            [dum, fr0] = max(abs(fft0));
            respamp = abs(fft0(fr0));                           %Calculate respiration amplitude
            respfreq = (fr0 - 1) / (length(WinResp) / sR);      %Calculate respiration frequency
            RVT_tsTTL(i) = respamp .* respfreq;                    %frequency times amplitude of respiration
        end
    end
    
    %Extrapolate beginning and end if NaNs
    if sum(~isnan(RVT_tsTTL))==0                %If there's no HR data at all
        RVT_tsTTL(:)=0;                         %Set to all zero
    end
    if isnan(RVT_tsTTL(1))                                          %If there's NaNs at the beginning
        stDatind = find(~isnan(RVT_tsTTL));                         %Get the first real datapoint
        RVT_tsTTL(1:stDatind(1)) = RVT_tsTTL(stDatind(1));          %Extrapolate to beginning
    end
    if isnan(RVT_tsTTL(end))                                        %If there's NaNs at the end
        endDatind = find(~isnan(RVT_tsTTL));                        %Get the last real datapoint
        RVT_tsTTL(endDatind(end):end) = RVT_tsTTL(endDatind(end)); %Extrapolate to end
    end
    
    %Collect the time shifted RVT vectors
    RVT_allTS(:,cTS_RVT) = RVT_tsTTL';
end %Loop over time shifts for RVT


%---------------------------------------------------------------
%Collect for output
%---------------------------------------------------------------

%Expand the cardiac and respiratory phase regressors to nth order
%Creating CPR and RPR output variables
F=HRPhase_TTL'*[1:fOrder];
CPR=[cos(F), sin(F)];                   %Cardiac phase regressors
F=[];                       %%added MvB
F=RespPhase_TTL'*[1:fOrder];              
RPR=[cos(F), sin(F)];                   %Respiratory phase regressors
%and combine all the HRF and RVT time shifts into NR
NR = [HRF_allTS,RVT_allTS];             %Nuisance regressors




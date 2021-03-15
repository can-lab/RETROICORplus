%--------------------------------------------------------------------------
%RETROICOR_defaults
%Default settings for RETROICORplus
%EJH 2011-16
%--------------------------------------------------------------------------

%Set Fourier order (3-5th order added to explained variance)
RETROICORplus_defaults.fOrder = 5;

%Window length (s) for heart rate calculations
RETROICORplus_defaults.HRFwinlen = 6;

%Time shifts (s) for heart rate (6 & 12 based on Shmueli et al. 2007, 10 based on van Buuren et al. 2009)
RETROICORplus_defaults.TS_HRF = [6,10,12];

%Window length (s) for RVT calculations
RETROICORplus_defaults.RVTwinlen = 9;

%Window length (s) for respiratory phase estimation
RETROICORplus_defaults.Respphasewinlen = 1;

%Time shifts for RVT (-1 & 5 based on Birn et al. 2006 and van Buuren et al. 2009)
RETROICORplus_defaults.TS_RVT = [-1,5];



%References:
%Shmueli, K., van Gelderen, P., de Zwart, J. A., Horovitz, S. G., Fukunaga, M., Jansma, J. M., & Duyn, J. H. (2007). Low-frequency fluctuations in the cardiac rate as a source of variance in the resting-state fMRI BOLD signal. NeuroImage, 38(2), 306?320.
%van Buuren, M., Gladwin, T. E., Zandbelt, B. B., van den Heuvel, M., Ramsey, N. F., Kahn, R. S., & Vink, M. (2009). Cardiorespiratory effects on default-mode network activity as measured with fMRI. Human Brain Mapping, 30(9), 3031?3042.
%Birn, R. M., Diamond, J. B., Smith, M. A., & Bandettini, P. A. (2006). Separating respiratory-variation-related fluctuations from neuronal-activity-related fluctuations in fMRI. NeuroImage, 31(4), 1536?1548.


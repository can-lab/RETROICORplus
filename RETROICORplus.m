%--------------------------------------------------------------------------
%RETROICORplus v2.3 (previously RETROICOR; v1.0-v1.1)
%
%RETROICORplus creates physiological noise regressors for fRMI data.
%
%usage:
%RETROICORplus(matFilename,RETROICORbegremscans,RETROICORendremscans, outputdir)
%
%matFilename: name of the matfile created by HeRa
%RETROICORbegremscans: how many EPIs are discarded at the beginning
%RETROICORendremscans:  how many EPIs are discarded at the end
%outputdir: folder where output regressors should be saved.
%
%This RETROICOR implementation is based on Glover et al's 2000 MRM paper.
%Initial versions of the implementation were created by Bas Neggers,
%Matthijs Vink, Thomas Gladwin, and Mariet van Buuren at Utrecht
%University. The current version was updated in collaboration with Mariet
%van Buuren.
%
%v2.2: Added a fix for a crash occurring when physiological recording is
%stopped only very briefly after scanning ends. This can cause crashes
%because the filtering will make the recording shorter, causing the TTL
%pulses to occur after the end of the recording.
%
%v2.3: Removed bug in cardiac phase regressor calculation
%
%EJH 2010-21
%--------------------------------------------------------------------------


function RETROICORplus(matFilename,RETROICORbegremscans,RETROICORendremscans, outputdir)

  
if ~nargin
    [matfile, matpath, irr] = uigetfile( ...
       {'*.mat','HeRa MAT files (*.mat)'}, ...
        'Pick a MAT file created by HeRa');
    matFilename = fullfile(matpath,matfile);
    
    outputdir = matpath;
    
    RETROICORbegremscans = inputdlg('Number of scans to remove at beginning','',1);
    RETROICORbegremscans = str2num(RETROICORbegremscans{1});

    RETROICORendremscans = inputdlg('Number of scans to remove at end','',1);
    RETROICORendremscans = str2num(RETROICORendremscans{1});

end

%--------------------------------------------------------------------------
%Get defaults
RETROICORplus_defaults_setup;

%--------------------------------------------------------------------------

%Load hera processed pulse data
heradata = load(matFilename);
%--------------------------------------------------------------------------

%Get sample rate
SR = heradata.matfile.settings.samplerate;

%Run rejection interpolation
heradata.matfile = ...
    RETROICORplus_interpolate_hera_reject(heradata.matfile,SR);

%Get all scan triggers, set to vertical vector
if size(heradata.matfile.markerlocs,2)>1
    scanTriggers = heradata.matfile.markerlocs';
else
    scanTriggers = heradata.matfile.markerlocs;
end

%Check if this is a continuous run
if sum(abs(diff(scanTriggers)-mean(diff(scanTriggers)))>2) >0 %TR is never more than 2 ms off mean
    error('This appears not to be a continuous recording')
end

%Run RETROICOR to create regressors
[CPR,RPR,NR]=RETROICORplus_calc(...
    scanTriggers,...                    %Scan trigger indices
    heradata.matfile.prepeaklocs,...    %Heart beat peak indices
    heradata.matfile.rawpulsedata,...   %Pulse data
    heradata.matfile.rawrespdata,...    %Respiration data
    SR,...                              %sample rate
    RETROICORplus_defaults);                %Defaults

%Save the result
R = [CPR,RPR,NR];
R = R(RETROICORbegremscans+1:end-RETROICORendremscans,:); %Remove omitted scans

[path,matname,EXT] = fileparts(matFilename);
outfile = fullfile(outputdir,[matname,'_RETROICORplus_regr.mat']);
save(outfile,'R');





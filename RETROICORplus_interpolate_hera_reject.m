%----------------------------------------------------
%RETROICOR_interpolate_hera_reject
%
%Interpolates heart rate data after editing in Hera
%to remove rejected periods by inserting new peaks.
%
%EJH2011-14
%----------------------------------------------------
function matfile=RETROICOR_interpolate_hera_reject(matfile,SR)
        
%Return if there's nothing to interpolate
if numel(matfile.prereject)==0
    return
end

%Get current pulse data and peak locations
cpulse=matfile.rawpulsedata;
cplocs=matfile.prepeaklocs;

%First create a rejection vector (because rejection periods can overlap)
recvec=zeros(size(matfile.rawpulsedata));
for crej=1:numel(matfile.prereject)
    %Get start/end sample of current rejection period
    crsts=matfile.prereject{crej}(1).*SR;
    %Set to one if lower than one
    crsts=max([crsts,1]);
    crens=matfile.prereject{crej}(2).*SR;
    warning off
    recvec(int32(crsts):int32(crens))=1;
    warning on
end

%Get starts and ends of rejection periods
rsts=find(diff([0,recvec,0])==1);
rens=find(diff([0,recvec,0])==-1);
if rens(end)>numel(recvec)
    rens(end)=numel(recvec);
end

%Loop over rejection periods to set rejections to NaN
for crej=1:numel(rsts)
    
    %Get last good peak before reject (if any)
    beforeRejPeaks = find(cplocs<rsts(crej)); %these are indices of cplocs
    if numel(beforeRejPeaks)> 0
        lastgood(crej) = beforeRejPeaks(end);
    else
        lastgood(crej) = 0;
    end
    
    %Get the first good peak after reject (if any)
    afterRejPeaks = find(cplocs>rens(crej)); %these are indices of cplocs
    if numel(afterRejPeaks)>0
        firstgood(crej)=afterRejPeaks(1);
    else %Set to number of elements in cplocs plus one, if there's none
        firstgood(crej)=numel(cplocs)+1;
    end
    
    %Remove all peaks in rejected period
    cplocs(lastgood(crej)+1:firstgood(crej)-1)=NaN;

end

allnewpeaklocs=[];

%Loop over rejection periods again to interpolate
for crej=1:numel(rsts)
    

    %Get an estimation on the mean IBI around that time
    
    %First get loc of middle of rejection period
    %cmloc=round(mean([rsts(crej),rens(crej)]));
    %Determine for every detected peak the distance to this
    %locdist=abs(cplocs-cmloc);
    %Take the closest 20 good beats, and calculate mean IBI
    %[irr,so]=sortrows(locdist');
    %meanIBIsamp=nanmean(diff(cplocs(min(so(1:20)):max(so(1:20)))));
    %Take the 20 closest IBIs to estimate a mean, first backward and then
    %forward
    
    %First go back 10 peaks, collect all non-rejected
    ccploc=lastgood(crej)-1; %start from one beat before last good
    iBack=0;
    cpreIBIsamp=[];
    while iBack<10 && ccploc>0 && ~isnan(cplocs(ccploc))
        %If nothing is rejected in this period
        if sum(recvec(cplocs(ccploc):cplocs(ccploc+1)))==0
            %Add the IBI (in samples)
            cpreIBIsamp=[cpreIBIsamp,cplocs(ccploc+1)-cplocs(ccploc)];
            iBack=iBack+1;

        end
        
        %go back one peak            
        ccploc=ccploc-1;
    
    end %while
    
    
    %Then go forward 10 peaks, collect all non-rejected
    ccploc=firstgood(crej); %get first good
    iForw=0;
    cpostIBIsamp=[];
    while iForw<10 && ccploc<numel(cplocs) && ~isnan(cplocs(ccploc+1))
        %If nothing is rejected in this period
        if sum(recvec(cplocs(ccploc):cplocs(ccploc+1)))==0
            %Add the IBI (in samples)
            cpostIBIsamp=[cpostIBIsamp,cplocs(ccploc+1)-cplocs(ccploc)];
            iForw=iForw+1;
        end
        
        %go forward one peak
        ccploc=ccploc+1;

    end %while
    
    %Calculate mean IBI in samples around rejection
    meanIBIsamp=mean([cpreIBIsamp,cpostIBIsamp]);
    
    %Now there a few possibilities
    %1 Ideal case: peaks around rejected period
    %2 No peaks before rejected period
    %3 No peaks after rejected period

    if lastgood(crej)>0 & firstgood(crej)<=numel(cplocs) %Case 1
        
        %Now take the time between the last good and first good
        %and determine how many peaks should be in there
        nsamphole = cplocs(firstgood(crej))-cplocs(lastgood(crej));
        nIBIshole = round(nsamphole/meanIBIsamp);
        newpeaklocs = cplocs(lastgood(crej)):(nsamphole/nIBIshole):cplocs(firstgood(crej));
        newpeaklocs = round(newpeaklocs(2:end-1));
        allnewpeaklocs = [allnewpeaklocs,newpeaklocs];
    
    elseif lastgood(crej)==0 & firstgood(crej)<=numel(cplocs) %Case 2
        
        %Now work backward from first good peak
        newpeaklocs = round([cplocs(firstgood(crej)):-meanIBIsamp:1]);
        newpeaklocs = newpeaklocs(2:end);
        allnewpeaklocs = [allnewpeaklocs,newpeaklocs];
        
    elseif lastgood(crej)>0 & firstgood(crej)>numel(cplocs) %Case 3
        
        %Now work forward from last good peak
        newpeaklocs = round([cplocs(lastgood(crej)):meanIBIsamp:numel(cpulse)]);
        newpeaklocs = newpeaklocs(2:end);
        allnewpeaklocs = [allnewpeaklocs,newpeaklocs];
    end
end

%Now remove NaNned peaklocs and add new ones and sort
newcplocs = sortrows([cplocs(~isnan(cplocs)),allnewpeaklocs]')';
newcptimes = newcplocs./SR;

%Calculate and store the pre-IBI timeseries
newibitimeseries = diff(newcptimes).*1000;

%And recalculate the IBI channel
newibichannel=[];
for cpeak = 1:length(newcplocs)
    %Create the pre-IBIchannel
    if cpeak>1
        newibichannel(newcplocs(cpeak-1):newcplocs(cpeak)) = ...
            newibitimeseries(cpeak-1);
    end
end

%Now fill up the start and end
if newcplocs(1)>1 
    newibichannel(1:newcplocs(1)-1) = newibitimeseries(1);
end

if newcplocs(length(newcplocs))<length(cpulse)
    newibichannel(newcplocs(length(newcplocs)):length(cpulse)) = ...
        newibitimeseries(length(newibitimeseries));
end

%Put new data back into structure
matfile.preibichannel = newibichannel;
matfile.prepeaklocs = newcplocs;
matfile.prepeakTimes = newcptimes;
matfile.preibitimeseries = newibitimeseries;


%Visualize
figure(1)
cla
hold on
for i=1:numel(newcplocs)
    cprepeak=newcplocs(i);
    plot([cprepeak,cprepeak],[0,4000],'b')
end
plot(cpulse,'r')
plot(recvec.*2000,'k')
plot(matfile.preibichannel,'r')
plot(newibichannel,'g')
hold off


        

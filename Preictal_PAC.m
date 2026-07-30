%%
clc;close all;
File = {'01_15'}; %names of files to run 
Seizure_start = [
1732]; % seizure start times of each file
Seizure_end = [
1772]; % seizure end times of each file
%seizureData = table(File,Seizure_start,Seizure_length);
numFiles = size(File,2);
medianPAC_allFiles = zeros(numFiles,2); 
medianPAC_maxChannelperFile = zeros(numFiles,2);
%%
for fileIndex = 1:numFiles
    if Seizure_start(fileIndex) > 600
% Load the dataset with changed channel label
load(['C:\Users\aditi\OneDrive\Documents\MATLAB\Examples\R2024b\Lab\chb',char(File(1,fileIndex)),'_reordered_resampled.mat']);  
if strcmp(File(fileIndex),'12_27') || strcmp(File(fileIndex),'12_28') || strcmp(File(fileIndex),'12_29')
    data = Rereference_EEG(reordered_record, reordered_hdr, 'bipolar');
    reordered_hdr.label = {'FP1F7' 'F7T7' 'T7P7' 'P7O1'	'FP1F3'	'F3C3' 'C3P3' 'P3O1' 'FP2F4' 'F4C4'	'C4P4'	'P4O2'	'FP2F8'	'F8T8'	'T8P8'	'P8O2'	'FZCZ'	'CZPZ'};
else
data = reordered_record;
end

%Check_EEG_Standardization(desired_channel_order, reordered_hdr);
fs = reordered_hdr.frequency(1); %save sampling rate
numChannels = size(data,1);
channel_order = reordered_hdr.label;

%% Find channel indices
nChan = length(reordered_hdr.label);
chanVec = nan(1, nChan);

for i = 1:nChan
    channelLoc = strcmp(reordered_hdr.label, channel_order{i});
    if all(~channelLoc)
        error('Channel %s not found in EEG header.', channel_order{i});
    end
    chanVec(i) = find(channelLoc, 1);
end
% Ensure EEG is channels x samples
if size(data,1) > size(data,2)
    warning(['EEG matrix appears to be samples x channels. ', ...
             'Transposing to channels x samples.']);
    data = data';
end
%% artifact detection
stdAbove  = 7.5;
buffer    = 0.9;
nArtChans = 1;
epochLength = 300;        % seconds
ARTepochsamples = 300*fs;
epochedData = data(:,(Seizure_start(fileIndex) - 600)*fs:Seizure_start(fileIndex,1)*fs);
nARTSamples = size(epochedData, 2);
nARTEpochs = ceil(nARTSamples / ARTepochsamples);
for e = 1:nARTEpochs
    % Sample range for this epoch
    startSample = (e-1)*ARTepochsamples + 1;
    endSample   = min(e*ARTepochsamples, nARTSamples);

    % Extract this epoch
    epochData = epochedData(:, startSample:endSample);
    % Here, we have to make a guardrail - if the size is small,
    % at the end of the signal - exclude them
    if e == nARTEpochs && size(epochData, 2) < 1*fs
        disp("The Artifact not excluded for the end segment.");
        continue;
    end
    [artifactInd, Artifacts] = Detect_Artifacts(epochedData, fs, stdAbove, buffer, nArtChans, numChans);
end

% find clean epochs
epochStart = Find_Clean_Indices(artifactInd, fs, epochLength);
epochStop  = epochStart + epochLength * fs - 1;
% nEpoch     = length(epochStart);
d = epochedData;
t = 1:fs:300*fs;
hold on
for i = 1:18
    if i < 9
        plot(t,d(:,i) + (9-i)*100)
    elseif i > 9
        plot(t,d(:,i) - (abs(9-i))*100)
    else
        plot(t, d(:,i))
    end
end
hold off
%% filter 
notchFilteredEEG = Filter_EEG(data, fs, 'notch');
amplitudeSignal  = Filter_EEG(notchFilteredEEG, fs, 'gamma_pac'); % 35–70 Hz
phaseSignal      = Filter_EEG(notchFilteredEEG, fs, 'beta'); 
%phaseSignal      = Filter_EEG(notchFilteredEEG, fs, 'delta'); 
broadband        = Filter_EEG(notchFilteredEEG, fs, 'broadband');
%% compute PAC 
% https://pmc.ncbi.nlm.nih.gov/articles/PMC9980882/#hbm26190-sec-0002 - direct modulation index
nEpoch = 1;
epochStop = Seizure_end(fileIndex,1) + 300;
epochStart = Seizure_start(fileIndex,1) - 300;
MI  = nan(nEpoch, nChan);
MVL = nan(nEpoch, nChan);

if epochStart < 0
    epochStart = 1;
end
amplitudeEpochEEG = amplitudeSignal(:,epochStart*fs:epochStop*fs);
phaseEpochEEG = phaseSignal(:,epochStart*fs:epochStop*fs);

subEpochLengthPAC = 5;
    [epochMI, epochMVL] = Calc_PAC(amplitudeEpochEEG, phaseEpochEEG, fs, subEpochLengthPAC);

    % Save median PAC value across sub-epochs
    MI(nEpoch, :)  = median(epochMI, 1);
    MVL(nEpoch, :) = median(epochMVL, 1);

% MI_median_all_channels = median(MI,2);
% MVL_median_all_channels = median(MVL,2);
 %store MI & MVL values for channel with highest value
% MI_median_max = max(MI);
% MVL_median_max = max(MVL);

%% PLV
%https://pmc.ncbi.nlm.nih.gov/articles/PMC3674231/#S6
[maxVal, index] = maxk(max(epochMI),2);
channel1 = index(1); %channels to apply PLV
channel2 = index(2);
subEpochLengthPLV = 5;
[nEpochs, startInd, stopInd] = Calc_Epoch_Indices(size(amplitudeEpochEEG,2),fs,subEpochLengthPLV);
PLV = nan(1,nEpochs);
%PLV
for subEpoch = 1:nEpochs
    startidx = startInd(subEpoch);
    endidx = stopInd(subEpoch);
    z1 = amplitudeEpochEEG(channel1,startidx:endidx) + j*hilbert(amplitudeEpochEEG(channel1,startidx:endidx));
    z2 = amplitudeEpochEEG(channel2,startidx:endidx) + j*hilbert(amplitudeEpochEEG(channel2,startidx:endidx));
    relative_phase = angle(z1.*(conj(z2))./(abs(z1).*abs(z2)));
    %relative_phase = angle(z1)-angle(z2);
    PLV(1,subEpoch) = abs(mean(exp(j*relative_phase))); %magnitude of expected value of relative phase (as complex number)
end
%{ 
t = 0:1/1000:15;
sin1 = sin(2*pi*t + pi/4*t);
sin2 = sin(2*pi*t + pi/2);
z1 = sin1 + j*hilbert(sin1);
z2 = sin2 + j*hilbert(sin2);
relative_phase = angle((z1.*(conj(z2)))./(abs(z1).*abs(z2)));
PLV = abs(mean(cos(relative_phase) + j*sin(relative_phase)));
%}
%% plot
%heatmap(MI_map)
%{
figure(fileIndex)
subplot(1,2,1)
    plot(MI,'*')
    title([fileNumber(4:5),'-',fileNumber(7:end),' MI'])
    xticks(1:18)
    xticklabels(channel_order)
    subplot(1,2,2)
    plot(MVL,'*')
    title([fileNumber(4:5),'-',fileNumber(7:end),' MVL'])
    xticks(1:18)
    xticklabels(channel_order)
%}
% medianPAC_allFiles(fileIndex,1) = MI_median_all_channels;
% medianPAC_allFiles(fileIndex,2) = MVL_median_all_channels;
% medianPAC_maxChannelperFile(fileIndex,1) = MI_median_max;
% medianPAC_maxChannelperFile(fileIndex,2) = MVL_median_max;

%plot
figure(fileIndex)
t = epochStart:subEpochLengthPAC:epochStop-subEpochLengthPAC;

subplot(3,1,1)
plot(t,median(epochMI,2))
ylabel('Modulation Index')
subplot(3,1,2)
plot(t,median(epochMVL,2))
ylabel('Mean Vector Length')
subplot(3,1,3)
t = epochStart:subEpochLengthPLV:epochStop-subEpochLengthPLV;
plot(t,PLV)
ylabel('Phase Locking Value')
for i = 1:3
    subplot(3,1,i)
    xlabel('Time (seconds)')
    xline(Seizure_start(fileIndex),'r')
    xline(Seizure_end(fileIndex),'r')
    if fileIndex == 2
    xline([1120 1438 1411 1745],'r')
    end
end
    end
end
%{
figure(fileIndex*2)
subplot(2,1,1)
t = ((0:size(epochMI,1)-1)*subEpochLengthPAC+(Seizure_start(fileIndex)-300+subEpochLengthPAC/2))*fs;
plot(t,median(epochMI,2))
%plot(t,max(epochMI,[],2))
xline(Seizure_start(fileIndex),'r')
xline(Seizure_end(fileIndex),'r')
%xline([1088 1411 1745 1120 1438 1764],'r')
subplot(2,1,2)
%plot(t,max(epochMVL,[],2))
plot(t,median(epochMVL,2))
%xline([1088 1411 1745 1120 1438 1764],'r')
xline(Seizure_start(fileIndex),'r')
xline(Seizure_end(fileIndex),'r')
%}

%% Plots
%{
% seizure length vs max pac per file
[sorted_length,index] = sort(Seizure_length);
figure(1)
title('Max PAC value vs seizure length')
plot(sorted_length,MIperFileMax(index))
xlabel('Seizure Length (sec)')
ylabel('Modulation Index (MI)')
figure(2)
plot(sorted_length,MVLperFileMax(index))
xlabel('Seizure Length (sec)')
ylabel('Mean Vector Length (MVL)')

% seizure length vs median pac per file
figure(3)
plot(sorted_length,MIperFile(index))
xlabel('Seizure Length (sec)')
ylabel('Modulation Index (MI)')
figure(4)
plot(sorted_length,MVLperFile(index))
xlabel('Seizure Length (sec)')
ylabel('Mean Vector Length (MVL)')
%}
% binned seizure length vs median MI

%{
numBins = 10;
figure(1)
[Seizure_length_no_outliers,TFrm,TFoutlier] = rmoutliers(Seizure_length);
MI_noOutliers = MIperFile;
MI_noOutliers(TFrm,:) = [];
[numSamplesPerBin,edges,seizureLengthBins] = histcounts(Seizure_length_no_outliers,numBins);

binnedMI = nan(max(numSamplesPerBin),numBins);
for binNum = 1:numBins
count = 1;
for i = 1:length(seizureLengthBins)
    if (seizureLengthBins(i,:) == binNum)
        binnedMI(count,binNum) = MI_noOutliers(i,1);
        count = count + 1;
    end
end
end

ttest = nan(numBins,2);
for binNum = 1:numBins-1
    [h,p] = ttest2(binnedMI(:,binNum),binnedMI(:,binNum+1));
    ttest(binNum,1) = h;
    ttest(binNum,2) = p;
end

%MI
boxchart(seizureLengthBins,MI_noOutliers);
sigstar([1,2],0.0019)
xticks(1:numBins);
xticklabels(compose('[%.0f, %.0f] ',edges(1:end-1)',edges(2:end)'))
xlabel('Seizure Length (sec)')
ylabel('Modulation Index (MI)')
for i = 1:10
text(i,4.5e-3,['n:',string(numSamplesPerBin(1,i))])
end

%MVL
boxchart(seizureLengthBins,MVL_noOutliers);
sigstar({[1,2],[6,7]},[0.0317,0.0108588679978469])
xticks(1:numBins);
xticklabels(compose('[%.0f, %.0f] ',edges(1:end-1)',edges(2:end)'))
xlabel('Seizure Length (sec)')
ylabel('Mean Vector Length (MVL)')
for i = 1:10
text(i,1,['n:',string(numSamplesPerBin(1,i))])
end

figure(2)
MVL_noOutliers = MVLperFile;
MVL_noOutliers(TFrm,:) = [];
[seizureLengthBins,edges] = discretize(Seizure_length_no_outliers,numBins);
boxchart(seizureLengthBins,MVL_noOutliers)
xticks(1:numBins);
xticklabels(compose('[%.0f, %.0f] ',edges(1:end-1)',edges(2:end)'))
xlabel('Seizure Length (sec)')
ylabel('Mean Vector Length (MVL)')

%}



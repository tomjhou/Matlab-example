% This script tests PLX import routines
% This script was downloaded from Plexon's website, then modified to make
% it more intuitive to new users.
disp('This sample script is provided for demonstration purposes. You should modify it to suit your own needs.');

% Open a plx file
% this will bring up the file-open dialog
StartingFileName = '';
%StartingFileName = 'C:\PlexonData\CM_Quickstart.plx';
[OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreThresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(StartingFileName);

disp([]);
disp(['Opened File Name: ' OpenedFileName]);
disp(['Version: ' num2str(Version)]);
disp(['Sampling frequency : ' num2str(Freq)]);
disp(['Comment : ' Comment]);
disp(['Date/Time : ' DateTime]);
disp(['Duration (seconds): ' num2str(Duration)]);
disp(['Num Pts Per Waveform : ' num2str(NPW)]);
disp(['Num Pts Pre-Threshold : ' num2str(PreThresh)]);
% some of the information is only filled if the plx file version is >102
if ( Version > 102 )
    if ( Trodalness < 2 )
        disp('Data type : Single Electrode');
    elseif ( Trodalness == 2 )
        disp('Data type : Stereotrode');
    elseif ( Trodalness == 4 )
        disp('Data type : Tetrode');
    else
        disp('Data type : Unknown');
    end
        
    disp(['Spike Peak Voltage (mV) : ' num2str(SpikePeakV)]);
    disp(['Spike A/D Resolution (bits) : ' num2str(SpikeADResBits)]);
    disp(['Slow A/D Peak Voltage (mV) : ' num2str(SlowPeakV)]);
    disp(['Slow A/D Resolution (bits) : ' num2str(SlowADResBits)]);
end   


% get counts of # of events for spike, waveform, event, and continuous data channels
[tscounts, wfcounts, evcounts, ad_counts] = plx_info(OpenedFileName,1);

% tscounts, wfcounts are indexed by (channel+1,unit+1)
% tscounts(:,ch+1) is the per-unit counts for channel ch
% sum( tscounts(:,ch+1) ) is the total wfs for channel ch (all units)
% [nunits, nchannels+1] = size( tscounts )
% To get number of nonzero units/channels, use nnz() function

% The next three variables help us look for first spike and continuous channels that have data (i.e. are
% non-empty).
firstSpikeChannel = -1;
firstSpikeChannelWithSortedUnits = -1;
firstSortUnit = -1; % This will be zero for unsorted spikes, 1 for sorted unit "a"
firstAD_Channel = -1;

disp('-');
% get spike channel names
[n, channames] = plx_chan_names(OpenedFileName);
disp(['Number of spike channels in this file : ' num2str(n)]);
for x = 1:n

    % plx_chan_names function returns null-padded string names. This can cause
    % problems later, e.g. when using names as graph titles.
    % The following line replaces all nulls with white space
    channames(channames == 0) = 32;
    
    % Display names of all channels that are non-empty.
    if tscounts(1,x+1) > 0
        % note the use of x+1, not x, in index for tscounts. Why did Plexon
        % have to make this so weird?
        disp(['Channel ' num2str(x) ' is named ' channames(x,:) ', and has ' num2str(tscounts(1,x+1)) ' spikes']);
    end

    if (firstSpikeChannel < 0) && (tscounts(1,x+1) > 0)
        % Find the first channel that has any spike data, whether sorted or
        % not
        firstSpikeChannel = x;
        firstSortUnit = 0;
    end
    if (firstSpikeChannelWithSortedUnits < 0) && (tscounts(2,x+1) > 0)
        % Find the first channel that has a SORTED unit "a"
        firstSpikeChannelWithSortedUnits = x;
        % Channel with sorted data will take precedence over channel with
        % only unsorted data
        firstSpikeChannel = x;
        firstSortUnit = 1;
    end
end

disp('-');
% get A/D channel names (EEG channels)
[n, ad_channames] = plx_adchan_names(OpenedFileName);
disp(['Number of continuous/EEG channels : ' num2str(n)]);
for x = 1:n

    % plx_adchan_names function returns null-padded string names. This can cause
    % problems later, e.g. when using names as graph titles.
    % The following line replaces all nulls with white space
    ad_channames(ad_channames == 0) = 32;

    % note that ad_counts is indexed normally, i.e. using x instead of x+1
    
    % Display names of all non-empty continuous channel
    if ad_counts(1,x) > 0
        disp(['EEG channel ' num2str(x) ' is named ' ad_channames(x,:) ', and has ' num2str(ad_counts(1,x)) ' points']);
    end

    % If this is the first non-empty continuous channel, then remember it.
    if (firstAD_Channel < 0) && (ad_counts(1,x) > 0)
            firstAD_Channel = x;
        else
    end
end

disp('-');

% Show spike timestamps
if firstSpikeChannel > 0

    % get some timestamps for channel 1 unit a
    [nts, ts] = plx_ts(OpenedFileName, 1, 1);
    
    disp(['Showing spike data from channel "' channames(firstSpikeChannel,:) '":']);

    % Display timestamps of first 10 samples
    numSamplesToShow = min([nts, 10]);
    disp(['First ' num2str(numSamplesToShow) ' samples have the following timestamps:']);
    for x = 1:numSamplesToShow
        disp(['  Sample ' num2str(x) ': ' num2str(ts(x)) ' sec']);
    end
    
    [nwf, npw, tswf, waves] = plx_waves(OpenedFileName, firstSpikeChannel, firstSortUnit);
    if nwf > 0

        if nwf >= 10
            % If there are >= 10 waves, then show first 10 waves.

            % Transpose to be compatible with what plot function expects in
            % 2D arrays.
            waves2 = waves';
            plot(waves2(1:NPW,1:10));
            % strtrim removes trailing blank space that plx_names appends
            title(['First 10 waves in channel "' strtrim(channames(firstSpikeChannel,:)) '"']);
        else
            % If there aren't >= 10 waves, then just plot first one
            plot( waves(1,1:NPW));
            % strtrim removes trailing blank space that plx_names appends
            title(['First waveform in channel "' strtrim(channames(firstSpikeChannel,:)) '"']);
        end
        
    end
else
    disp('No spike channels with data were found in this file.');
end
    
% Show EEG/continuous data
if firstAD_Channel > 0
    
    % get all data from specified EEG/continuous channel
    [adfreq, nad, tsad, fnad, ad] = plx_ad(OpenedFileName, firstAD_Channel - 1);
    if nad > 0
        % get just a subset of a/d data for sample plot
        [adfreq, nadspan, adspan] = plx_ad_span(OpenedFileName, firstAD_Channel - 1, 1,1000);

        disp(['Showing continuous data from channel "' ad_channames(firstAD_Channel,:) '":']);


        % Send plot to a new window
        figure;
        % Plot first 1000 samples in red
        plot( adspan, 'r' );
        title(['First ' num2str(nadspan) ' samples of continuous/EEG channel "' strtrim(ad_channames(firstAD_Channel,:)) '"'], 'interpreter', 'none');
    end

end

% Close all open files
plx_close('');

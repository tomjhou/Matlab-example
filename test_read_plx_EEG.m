% This script tests PLX import routines
%
% This script was downloaded from Plexon's website, then modified to make
% it more intuitive to new users.
disp(' ');

% Open a plx file
[OpenedFileName, Version, SamplingFreq, Comment, Trodalness, PointsPerWaveform, PointsPerWaveformPreThreshold, ...
    SpikeChannelMaxRangeMillivolts, SpikeChannelResolutionBits, ContinuousChannelMaxRangeMillivolts, ContinuousChannelResolutionBits, ...
    DurationSeconds, DateTimeString] = plx_information('');

disp(' ');
disp(['File Name: ' OpenedFileName]);
disp(['Version: ' num2str(Version)]);
disp(['Sampling frequency : ' num2str(SamplingFreq)]);
disp(['Comment : ' Comment]);
disp(['Date/Time : ' DateTimeString]);
disp(['File Duration (seconds): ' num2str(DurationSeconds)]);

disp(['Continuous Channel Voltage Range (positive only) (mV) : ' num2str(ContinuousChannelMaxRangeMillivolts)]);
disp(['Continuous channel A/D Resolution (bits) : ' num2str(ContinuousChannelResolutionBits)]);


% get counts of # of events for spike, waveform, event, and continuous data channels
[spike_counts, waveform_counts, event_counts, continuous_sample_counts] = plx_info(OpenedFileName,1);


firstContinuous_Channel_Wire = -1;
firstContinuous_Channel_Index = -1;



disp(' ');
% Get EEG channel wire #. Normally the wires are sequential (1, 2, 3 ...)
% but there could be non-consecutive channels, or map could starts at 0.
% Wire number is used to retrieve raw data.
[n, continuous_channel_wire] = plx_ad_chanmap(OpenedFileName);
% get EEG channel names
[n, continuous_channames] = plx_adchan_names(OpenedFileName);
fprintf('Number of continuous/EEG channels in file header (not all may have data) : %d\n', n);
for x = 1:n

    % plx_adchan_names function returns null-padded string names.
    % We replace nulls with white space to avoid problems later, e.g. when
    % names are used as graph axis/legend/title text.
    continuous_channames(continuous_channames == 0) = 32;

    % Display names of all non-empty continuous channel
    if continuous_sample_counts(x) > 0
        fprintf('EEG channel %d is named %s, and has %d points\n', x, continuous_channames(x,:), continuous_sample_counts(x));
    end

    % If this is the first non-empty continuous channel, then remember it.
    if (firstContinuous_Channel_Wire < 0) && (continuous_sample_counts(x) > 0)
            firstContinuous_Channel_Index = x;
            firstContinuous_Channel_Wire = continuous_channel_wire(x);
        else
    end
end

    
% Show EEG/continuous data
if firstContinuous_Channel_Index > 0
    
    % get all data from specified EEG/continuous channel. Note that the
    % channel index here is ZERO-based, so we have to subtract one from the
    % ONE-based channel number in firstAD_Channel.

    fprintf('\nSelect option for plotting EEG/continuous data:\n');
    fprintf('1: Read and plot all data from first EEG/continuous channel (not recommended for very large files)\n');
    fprintf('2: Read and plot first 1000 samples only\n');
    answer = input('Please select option: ');

    % Get channel gains. 
    [n,continuous_channel_gains] = plx_adchan_gains(OpenedFileName);
    
    % Calculate scaling factors for converting raw values to voltages
    % The following is the # of samples between lowest and highest voltage.
    % Typically 2^16 = 65536.
    sample_range = 2 ^ ContinuousChannelResolutionBits;
    % The following is the difference between highest and lowest possible
    % voltage. Note the multiplication by 2, because this range is
    % two-sided, i.e. includes negative numbers as well.
    voltage_range = 2 * ContinuousChannelMaxRangeMillivolts / 1000 / continuous_channel_gains(firstContinuous_Channel_Index);

    if answer ~= 1 && answer ~= 2
        fprintf('Invalid option\n');
        plx_close('');
        return;
    end
    
    if answer == 1
        
        % Read all data. If file is very large, this could suck up a lot of
        % memory.
    
        % plx_ad_span retrieves all data into one large vector. If there are
        % any acquisition gaps, they will be padded with zeros.
        [continuous_sample_freq, total_num_datapoints, raw_ad_values] = plx_ad_span(OpenedFileName, firstContinuous_Channel_Wire, 1, continuous_sample_counts(firstContinuous_Channel_Index));

        % One disadvantage of plx_ad_span is that you have to already know the
        % final timestamp in order to get all data (it is the last argument passed).
        % If you instead use plx_ad(), you don't have to know the data length. But it
        % may retrieve data in several chunks, which you'll have to concatenate
        % back together. The number of chunks is somewhat
        % machine-dependent, and if there are gaps, you'll have to pad them
        % yourself. Your choice.
    %    [adfreq, total_num_datapoints, timestamps, fragment_lengths, raw_ad_values] = plx_ad(OpenedFileName, firstAD_Channel - 1);

    elseif answer == 2
        
        % get first 1000 samples only.
        [continuous_sample_freq, total_num_datapoints, raw_ad_values] = plx_ad_span(OpenedFileName, firstContinuous_Channel_Wire, 1,1000);

    else
        
        fprintf('Invalid option\n');
        
    end
    
    fprintf('Showing continuous data from channel "%s"\n', continuous_channames(firstContinuous_Channel_Index,:));

    % Convert raw values to voltages
    continuous_voltages = raw_ad_values / sample_range * voltage_range;

    % Send plot to a new window
    figure;

    % Calculate timestamps in seconds
    time_values = (1:length(continuous_voltages)) / continuous_sample_freq;

    % Plot first 1000 samples as voltage versus time (seconds) in red
    plot(time_values, continuous_voltages, 'r' );
    title(['Showing ' num2str(total_num_datapoints) ' samples of continuous/EEG channel "' strtrim(continuous_channames(firstContinuous_Channel_Index,:)) '"'], 'interpreter', 'none');
    xlabel('Time (seconds)');
    ylabel('Voltage');

end


% Close all open files
plx_close('');

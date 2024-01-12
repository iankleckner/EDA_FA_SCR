function SCR_Results = adaptiveEDA(X_sec, Y_uS, RAP_threshold, use_debug_plots)  

% adaptiveEDA - apply fixed + adaptive (FA) thresholding to electrodermal 
%   activity analysis (EDA) of skin conductance responses (SCR) 
%   (see Kleckner et al. *). Returns several metrics of SCRs using the FA 
%   thresholding method
% 
%   Requires "Signal Processing Toolbox" for function "findpeaks"
%
%   Input:
%       1) X - time vector of EDA signal in sec
%       2) Y - EDA signal vector in uS
%       3) RAP_threshold - Selected response amplitude percent (RAP)
%       threshold selected in percent (ex. 5 = 5%)
%       4) OPTIONAL use_debug_plots - true/false or 1/0 to create a figure
%       showing the EDA data and each SCR trough, peak, and HR point
%
%   Output:
%       1) SCR onset times (s)
%       2) SCR peak times (s)
%       3) SCR peak values (uS)
%       4) SCR amplitude value (uS)
%       5) SCR half recovery point times (s)
%       6) SCR half recovery point value (uS)
%       7) SCR total number 
%       8) Average SCL (uS)
%      
%   Output Notes: 
%       Half-recovery point is defined as the first point below the half
%       recovery threshold between the peak of the SCR and the peak of the 
%       next SCR. If that point does not exist, an NaN is saved in the 
%       place of that time index.
%
%       SCL/Skin conductance level is defined as the average of EDA 
%       signal provided as the input to this function, excluding segments 
%       of data that contain SCRs (either from trough to half recovery 
%       point or trough to peak if no half recovery point exists). 
%
%   Maya Barton-Zuckerman
%   PhD Student 
%   Northeastern University 
%   barton-zuckerman.m@northeastern.edu
%
%   Ian Kleckner, PhD, MPH
%   Associate Professor
%   University of Maryland Baltimore
%   ian.kleckner@umaryland.edu
%
%------CHANGELOG-----------------------------------------------------------
%   2023/12/01  MBZ wrote first draft with IK
%   2024/01/12  IK updated to save as results structure

% User options/settings: 
%   Recommended parameters:
%       peak_prominence_uS_MINIMUM        = 0.01
%       min_time_from_peak_to_trough_sec  = 1
%       max_time_from_peak_to_trough_sec  = 3
%
%   debug_plots - creates plot of input data and resulting EDA annotations 
%   (peak, trough, and half recovery points) for visual inspection
% -------------------------------------------------------------------------
peak_prominence_uS_MINIMUM        = 0.01; 
min_time_from_peak_to_trough_sec  = 1; 
max_time_from_peak_to_trough_sec  = 3; 

if( nargin < 4 )
    use_debug_plots = false; 
end
% -------------------------------------------------------------------------

% Calculate sampling period based on time vector provided
sampling_period                   = X_sec(2)-X_sec(1);

% Detect the peaks in the EDA signal
%   peak_EDA_uS         Local maxima, returned as a vector of signal values.
%   peak_time_sec       Peak locations, returned as a vector
[peak_EDA_uS_FP, peak_time_sec_FP, ~, ~] = ...
    findpeaks(Y_uS, X_sec, 'MinPeakProminence', peak_prominence_uS_MINIMUM);

% Detect troughs in the EDA signal based on peak detection
%   trough_struct.INPUT_PEAK_HAS_VALID_RISE_TIME   Array w/ 0 or 1 representing valid peaks/troughs
%   trough_struct.EDA_uS_x_Peaknum                 8 sec of EDA signal before each peak
%   trough_struct.EDA_peak_array                   SCL value at the peak 
%   trough_struct.EDA_trough_array                 SCL value at the trough
%   trough_struct.rise_time_sec_array              Time between peak/trough in sec
%   trough_struct.rise_time_locs                   Indices of valid peaks/troughs
%   trough_struct.raw_peak_latency                 Time between peak/trough in samples  
trough_struct = EDA_trough_det( 1/sampling_period, peak_time_sec_FP, X_sec, Y_uS, ...
    min_time_from_peak_to_trough_sec, max_time_from_peak_to_trough_sec); 

% Select the proper peaks and get relevant info about peaks/troughs, etc.
% for later calculations for SCR metrics
peak_EDA_uS         = peak_EDA_uS_FP( trough_struct.INPUT_PEAK_HAS_VALID_RISE_TIME );
trough_EDA_uS       = trough_struct.EDA_trough_array( trough_struct.INPUT_PEAK_HAS_VALID_RISE_TIME )';
peak_time_sec       = peak_time_sec_FP( trough_struct.INPUT_PEAK_HAS_VALID_RISE_TIME )';
rise_time_sec       = trough_struct.rise_time_sec_array( trough_struct.INPUT_PEAK_HAS_VALID_RISE_TIME);

% Calculation of prominence using PRECEDING trough
peak_prominence_uS  = peak_EDA_uS - trough_EDA_uS;

% Now remove any peaks that have a prominence of < threshold,
PEAK_IS_VALID = peak_prominence_uS > peak_prominence_uS_MINIMUM;
peak_EDA_uS         = peak_EDA_uS(PEAK_IS_VALID);
trough_EDA_uS       = trough_EDA_uS(PEAK_IS_VALID);
peak_time_sec       = peak_time_sec(PEAK_IS_VALID);
rise_time_sec       = rise_time_sec(PEAK_IS_VALID);
peak_prominence_uS  = peak_prominence_uS(PEAK_IS_VALID);

% Calculation of response amplitude percent threshold
peak_EDA_prominence_percent = 100*peak_prominence_uS ./ trough_EDA_uS;

% Now remove any peaks that have response amplitude percent threshold < threshold
PEAK_MEETS_THRESHOLD = peak_EDA_prominence_percent > RAP_threshold;
peak_EDA_uS         = peak_EDA_uS(PEAK_MEETS_THRESHOLD);
peak_time_sec       = peak_time_sec(PEAK_MEETS_THRESHOLD);
rise_time_sec       = rise_time_sec(PEAK_MEETS_THRESHOLD);
peak_prominence_uS  = peak_prominence_uS(PEAK_MEETS_THRESHOLD);
trough_EDA_uS       = trough_EDA_uS(PEAK_MEETS_THRESHOLD);

% Calculate onset time of SCRs
trough_time_sec     = peak_time_sec - rise_time_sec;

% Calculate half recovery time (time from peak to time of return to 50% of peak prominence) 
peak_half_recovery_thresh = peak_EDA_uS - (peak_prominence_uS ./ 2); 
peak_half_recovery_sec = zeros(size(peak_time_sec));
peak_half_recovery_uS = zeros(size(peak_time_sec));
for peak = 1:length(peak_EDA_uS)
    peak_time_val = peak_time_sec(peak);
    peak_ind = find(X_sec==peak_time_val);
    if peak < length(peak_EDA_uS)
        next_peak_time_val = peak_time_sec(peak + 1);
        next_peak_ind = find(X_sec==next_peak_time_val);
        search_window = Y_uS(peak_ind:next_peak_ind);
    else 
        search_window = Y_uS(peak_ind:end);
    end
    half_recovery_window_ind = find(search_window <= peak_half_recovery_thresh(peak), 1, 'first');
    if isempty(half_recovery_window_ind)
        peak_half_recovery_sec(peak) = NaN; 
    else
        half_recovery_ind = peak_ind + half_recovery_window_ind;
        half_recovery_sec = X_sec(half_recovery_ind);
        peak_half_recovery_sec(peak) = half_recovery_sec;
        peak_half_recovery_uS(peak) = search_window(half_recovery_window_ind);
    end 
end

% Calculate SCL
data_for_SCL = Y_uS;
for SCR = 1:length(peak_time_sec)
    start_ind = find(X_sec == trough_time_sec(SCR));
    if isnan(peak_half_recovery_sec(SCR))
        end_ind = find(X_sec == peak_time_sec(SCR));
    else
        end_ind = find(X_sec == peak_half_recovery_sec(SCR));
    end
    data_for_SCL(start_ind:end_ind) = NaN;
end

% Set output variables
SCR_onset_time          = trough_time_sec;
SCR_onset_value         = trough_EDA_uS;
SCR_peak_time           = peak_time_sec;
SCR_peak_value          = peak_EDA_uS;
SCR_amplitude_value     = peak_prominence_uS;
SCR_half_recovery_time  = peak_half_recovery_sec;
SCR_half_recovery_value = peak_half_recovery_uS;
SCR_total_num           = length(SCR_peak_time);
SCL_avg                 = mean(data_for_SCL,"omitnan");

% Save results to output structure
SCR_Results.SCR_onset_time              = SCR_onset_time;
SCR_Results.SCR_onset_value             = SCR_onset_value;
SCR_Results.SCR_peak_time               = SCR_peak_time;
SCR_Results.SCR_peak_value              = SCR_peak_value;
SCR_Results.SCR_amplitude_value         = SCR_amplitude_value;
SCR_Results.SCR_half_recovery_time      = SCR_half_recovery_time;
SCR_Results.SCR_half_recovery_value     = SCR_half_recovery_value;
SCR_Results.SCR_total_num               = SCR_total_num;
SCR_Results.SCL_avg                     = SCL_avg;

% Create debugging plots
if (use_debug_plots)
    plot(X_sec,Y_uS,'-k') % plot raw
    hold on 
    scatter(SCR_peak_time,SCR_peak_value,'*r'); % plot SCR peaks
    hold on
    scatter(SCR_onset_time,SCR_onset_value,'og') % plot SCR trough
    hold on
    scatter(SCR_half_recovery_time,SCR_half_recovery_value,'om') % plot half recovery point

    % Set plot titles, axis labels, and legend 
    title('Fixed-Adaptive Thresholding EDA analysis');
    xlabel('Time (s)');
    ylabel('Electrodermal Activity (uS)');
    legend({'EDA signal','SCR Peaks','SCR Troughs','SCR Half-Recovery'},'Location','northwest');
end

end


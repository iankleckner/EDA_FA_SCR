function [trough_struct]=EDA_trough_det(hz,peak_time_sec,X,Y,min_t,max_t)

% James Heathers
% Some comments and variable name changes by Ian Kleckner
% 2019/07/22
%
% INPUT
% hz is the sampling rate
% peak_time_sec is X values of the detected peaks
% X is the time vector
% Y is the SCR values, corrected. May included NaNs.
% min_t = 1 sec recommended min value for SCR detection
% max_t = 3 sec recommended max value for SCR detection 

% 8 seconds, time for establishing SCR width taken pre-peak - unlikely to be
% longer than that.

% Lazy tag!
%#ok<*AGROW>

        time_delay=8*hz; 
        data_EDA_x_peaknum = [];
        distance_peak_to_trough_samples = [];
        
        % Iterate through peaks peak = 1,2,3,...
        for peak=1:length(peak_time_sec)
            
            % Find time of peak and its index in the time array
            time_val=peak_time_sec(peak);
            time_loc=find(X==time_val);
            
            % Save the EDA data stream leading up to the peak
            if( time_loc-time_delay > 0 )
                if sum(isnan(Y(time_loc-time_delay:time_loc)))==0
                    data_EDA_x_peaknum(:,peak)=Y(time_loc-time_delay:time_loc); 
                else
                    data_EDA_x_peaknum(:,peak)=zeros(time_delay+1,1);
                end
            else
                data_EDA_x_peaknum(:,peak)=zeros(time_delay+1,1);
            end
        end
        
        
        if( ~isempty( data_EDA_x_peaknum ) )
            
            % For each EDA segment leading up to each peak, segment=1,2,3,...
            Nsegments = size(data_EDA_x_peaknum,2);
            for segment=1:Nsegments
                % If there are actually data here preceding the peak...
                %  pp is a matrix: EDA_uS x PEAKNUM (1,2,3...)

                if data_EDA_x_peaknum(1,segment) > 0

                    % Reverse the order of elements so the past is toward the
                    % right (larger indices)
                    single_SCR = flipud(data_EDA_x_peaknum(:,segment));

                    % Look for trough (2nd derivative > 0)
                    if min( find(diff(sign(diff(single_SCR))))+1 ) > 0
                        % Time from peak to trough in samples
                        distance_peak_to_trough_samples(segment)=min(find(diff(sign(diff(single_SCR))))+1); 

                        % SC value at trough
                        EDA_trough_array(segment) = single_SCR( distance_peak_to_trough_samples(segment) );

                        % SC value at peak
                        EDA_peak_array(segment) = single_SCR(1);
                    end
                end
            end
        end
        
        if( ~isempty( distance_peak_to_trough_samples ) )
            % Is the rise time within the desired time frame of the peak (e.g.,
            % within 1-3 sec)
            [~,good_rise_times] = find(distance_peak_to_trough_samples >= min_t*hz &  distance_peak_to_trough_samples <= max_t*hz);

            trough_struct.INPUT_PEAK_HAS_VALID_RISE_TIME = distance_peak_to_trough_samples >= min_t*hz &  distance_peak_to_trough_samples <= max_t*hz;


            trough_struct.EDA_uS_x_Peaknum      = data_EDA_x_peaknum; % bundle of SCR curves
            trough_struct.EDA_peak_array        = EDA_peak_array;
            trough_struct.EDA_trough_array      = EDA_trough_array;
            trough_struct.rise_time_sec_array   = distance_peak_to_trough_samples / hz;


            trough_struct.rise_time_locs        = good_rise_times; % acceptable rise times
            trough_struct.raw_peak_latency      = distance_peak_to_trough_samples; % times taken to reach SCR trough

        else
            % No peaks and troughs to return info on
            trough_struct.INPUT_PEAK_HAS_VALID_RISE_TIME = [];
            trough_struct.EDA_uS_x_Peaknum      = [];
            trough_struct.EDA_peak_array        = [];
            trough_struct.EDA_trough_array      = [];
            trough_struct.rise_time_sec_array   = [];
            trough_struct.rise_time_locs        = [];
            trough_struct.raw_peak_latency      = [];
        end

end

% save output, not used
% nsscr_filenm=[filename_Q_short(1:12) '_' num2str(p) '.mat'];
% save(nsscr_filenm,'trough_struct');
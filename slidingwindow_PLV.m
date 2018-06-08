% This code calculates PLV according to this approach:
%
% Lachaux, Jean?Philippe, Eugenio Rodriguez, Jacques Martinerie, and Francisco J Varela. 
% ?Measuring Phase Synchrony in Brain Signals.? Human Brain Mapping 8, no. 4 (January 1, 1999): 194?208. 
%
% with the additional change that it is performed in a sliding window.
% Note that phase lag consistency is defined across the time point within
% each window rather than as traditionally defined across trials. 
%
% Requires FieldTrip for band-pass filtering
%
% INPUTS:
%   	rootdir = string of path
%       subjects_indx = list of subect numbers
%
% Sepideh Sadaghiani 2018

function outfilepath = slidingwindow_PLV(rootdir, subjects_indx)

%- for band filter
all_bandstring{1} = '[1 3]';
all_bandstring{2} = '[4 7]';
all_bandstring{3} = '[8 10]';
all_bandstring{4} = '[10 12]';
all_bandstring{5} = '[13 25]';
all_bandstring{6} = '[26 45]';

%- for file naming
all_out_prefix{1} = 'delta';
all_out_prefix{2} = 'theta';
all_out_prefix{3} = 'lowalpha';
all_out_prefix{4} = 'highalpha';
all_out_prefix{5} = 'beta';
all_out_prefix{6} = 'gamma';

%% ================= calculate PLV subject-by-subject =====================
for subj = 1:length(subjects_indx)

    %% specify preprocessed input data
    clear compclean_data plv
     
    preprocfile = fullfile(rootdir, 'DATA', ['subject' num2str(subj)], 'functional/EEG', ['Preproc_subject' num2str(subj) '.mat'])
    
    if exist(preprocfile, 'file');
        
        disp(['subject' num2str(subj)]);
        
        load(preprocfile)                       % load preprocessed data and chnInfo
        chnData = compclean_data.trial{1};    	% actual data time courses (numchannel x time points)
        
        numChannels = size(chnData,1);
        
        for thisBand = 1:length(all_bandstring)       

            %% Band-pass filtering data to desired oscillation band
            bandstring  = all_bandstring{thisBand};
            disp(['band-pass filtering for band ' bandstring]);
            srate = compclean_data.fsample;
            freqband = eval(bandstring);                    % band of interest in Hz
            FIRorder = round(4*srate/median(freqband));     % order of the FIR filter. A useful rule of thumb can be to include about 4 to 5 cycles of the desired signal. For example,
                                                            % FIRorder = 100 for eeg data sampled at 250 Hz (=4ms cycle) corresponds to 400 ms and contains ~4 cycles of alpha band (10 Hz).
            filteredData = ft_preproc_bandpassfilter(chnData, srate, freqband, FIRorder, 'fir');
            
            %% apply hilbert transform to all channels. output is complex numbers (numchannel x time points)
            disp('hilbert transforming');
            transformedData = transpose(hilbert(transpose(filteredData)));  % identical to: ft_preproc_hilbert(dat,'no') in FieldTrip
            envelope_hilbert= abs(transformedData);
            transformedData = angle(transformedData);                       % into radiants
            
            %% calculate PLV
            disp('calculating PLV');
            % specify sliding window
            windowSize = 10*srate;      % timepoints: s*Hz
            timeStep   = 2*srate;       % time resolution in s*Hz
            allWindowCentre = windowSize/2 : timeStep : size(transformedData,2)-windowSize/2;
            plv = zeros(numChannels, numChannels, length(allWindowCentre));
            
            for channelCount = 1:numChannels-1
                
                for compareChannelCount = channelCount+1:numChannels
                    
                    differenceChannelData = transformedData(channelCount,:) - transformedData(compareChannelCount,:);
                    
                    for windowCount = 1:length(allWindowCentre)
                        thisWindow = allWindowCentre(windowCount)-(windowSize/2-1) : allWindowCentre(windowCount)+(windowSize/2);
                        
                        plv(channelCount, compareChannelCount, windowCount) = ...
                            abs(sum(exp(1i*(differenceChannelData(thisWindow))))) / length(thisWindow);
                    end
                end
            end
            
            %% global average of all pairwise PLVs
            disp('averaging across channel pairs for global PLV');
            globalplv = zeros(1,length(allWindowCentre));
            for windowCount = 1:length(allWindowCentre)
                thisWindowPlv = reshape(plv(:,:,windowCount), 1,numChannels^2);
                globalplv(windowCount) = mean(thisWindowPlv(find(thisWindowPlv)));
            end
            
            %% save data for this subject and this session
            outdir = fullfile(rootdir, 'DATA', ['subject' num2str(subj)], 'functional/EEG', 'PLV_values');
            if ~exist(outdir, 'dir'), mkdir(outdir), end
            outfilepath = fullfile(outdir, ['PLV_' all_out_prefix{thisBand} '_subject' num2str(subj) '.mat']);
            save(outfilepath, 'allWindowCentre', 'plv', 'globalplv', 'transformedData', '-v7.3', 'envelope_hilbert');
        end
       
    end
    
end


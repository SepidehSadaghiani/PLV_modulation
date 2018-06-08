% Simply loads individual subjects global PLV time course into an allsj
% cell, then normalizes the PLV time course of each band by that of all
% bands, and finally saves the mean, std and range of the PLV time course
% INPUTS:
%       rootdir = string of path
%       subjects_indx = list of subect numbers
%       do_norm = boolean. 1 means devide by all bands including band at hand.
%
% Sepideh Sadaghiani 2018

function allsubject_savepath = modulation_globalPLV(rootdir, subjects_indx, do_norm)

%% bands
all_out_prefix{1} = 'delta';
all_out_prefix{2} = 'theta';
all_out_prefix{3} = 'lowalpha';
all_out_prefix{4} = 'highalpha';
all_out_prefix{5} = 'beta';
all_out_prefix{6} = 'gamma';

%% load all sjs data into one large cell
globalPLV_allbands = cell(max(subjects_indx), 1);

for thisSubj = subjects_indx
    thissj_globalPLV_allbands = [];
    for this_band = 1:length(all_out_prefix)
        clear globalplv
        %- specify subject's input file
        thissj_PLV_file = fullfile(rootdir, 'DATA', ['subject' num2str(thisSubj)], 'functional', 'EEG', 'PLV_values_recalculated',...
            ['PLV_' all_out_prefix{this_band} '_subject' num2str(thisSubj) '_recalclated_scriptMay6copy.mat']);
        load(thissj_PLV_file)       
        thissj_globalPLV_allbands(this_band,:) = globalplv;        
    end  
    globalPLV_allbands{thisSubj} = thissj_globalPLV_allbands;
end

%% save mean, std and range after normalizing PLV by PLV of all band for each time window
mean_globalPLV_allsubj	= zeros(max(subjects_indx),1);
std_globalPLV_allsubj	= zeros(max(subjects_indx),1);
range_globalPLV_allsubj	= zeros(max(subjects_indx),1);

for this_band = 1:length(all_out_prefix)
    %- specify output for allsj file for this band
    out_folder = fullfile(rootdir, 'allsj_EEG_results', 'PLV_summary_values');
    if do_norm
        outfile = ['PLV_' all_out_prefix{this_band} '_allsj_norm_globalPLV_mean_std.mat'];
    else
        outfile = ['PLV_' all_out_prefix{this_band} '_allsj_nonnorm_globalPLV_mean_std.mat'];
    end
    allsubject_savepath = fullfile(out_folder, outfile);
    
    %- if requested, normalize. Then calculate mean, std, range
    for thisSubj = subjects_indx
        if do_norm
            tmp_PLV = globalPLV_allbands{thisSubj}(this_band,:)./mean(globalPLV_allbands{thisSubj});
        else
            tmp_PLV = globalPLV_allbands{thisSubj}(this_band,:);
        end
        mean_globalPLV_allsubj(thisSubj)	= mean(tmp_PLV);
        std_globalPLV_allsubj(thisSubj)     = std(tmp_PLV);
        range_globalPLV_allsubj(thisSubj)   = max(tmp_PLV)-min(tmp_PLV);
    end
    
    % save allsj output
    save(allsubject_savepath, 'mean_globalPLV_allsubj', 'std_globalPLV_allsubj', 'range_globalPLV_allsubj');
end


function [] = read_persist_logfiles(subject,data_dir)
% Jamil P. Bhanji 
% bhanji@utexas.edu
% 12/18/2012

%NOTE: some of the plots require functions written by other Matlab users:
%   add_errorbars.m: written by Petr Janata
%   xticklabel_rotate.m: written by Brian FG Katz, based on
%   xticklabel_rotate90.m (written by Denis Gilbert)
%   xticklabel_rotate is available on Matlab central file exchange

% 1. read in one e-prime text logfile for one subject. 
% 2. output fsl onset files for each model for each run
%%%%% not working yet 3. Append a line to output file w/ summary of behavioral data
% call with subject as 3 digit string, e.g. '001', and data_dir as folder
% with the e-prime txt file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EXP_DIR = '/mnt/delgadolab/jamil/persist/';

%if length(subject==4) subject = subject(2:end); end

WRITE_BEHAV = 1;

subject_dir = fullfile(EXP_DIR,sprintf('sub%s',subject));
if ~exist('data_dir','var') data_dir=fullfile(subject_dir,'behav'); end
DATADIR = data_dir;  
DATAFILE_PRE = 'PersTasknewcomp_fmri_1234-';
DATAFILE_SUF = '-1.txt';



num_runs = 4;
files_to_read = 1;
trials_per_run = 88;

%vars and masks for behavioral data (length = 88 events per round), for use w/
persist_data = zeros(trials_per_run,num_runs);
switch_data = zeros(trials_per_run,num_runs);
rt_data = zeros(trials_per_run,num_runs);

%respcond_mask = false(files_to_read * trials_per_run, 6); %columns correspond to 6 event types that require a response
% cond001 first choice of a round
% cond002 choice after uncontrollable failure
% cond003 choice after controllable failure
% cond004 positive cue
% cond005 uncontrollable cue
% cond006 controllable cue
choice_onset_list = zeros(trials_per_run, num_runs); 
cue_onset_list = zeros(trials_per_run, num_runs);
fbk_onset_list = zeros(trials_per_run, num_runs);
goal_onset_list = zeros(trials_per_run, num_runs);
choice_duration_list = zeros(trials_per_run, num_runs); choice_duration_list(:,:) = 2;
cue_duration_list = zeros(trials_per_run, num_runs); cue_duration_list(:,:) = 2;
fbk_duration_list = zeros(trials_per_run, num_runs); fbk_duration_list(:,:) = 2;
goal_duration_list = zeros(trials_per_run, num_runs); goal_duration_list(:,:) = 2;
cuefbk_duration_list = zeros(trials_per_run, num_runs); cuefbk_duration_list(:,:) = 4;
weights_list = ones(trials_per_run, num_runs);

firstchoice_onset_mask = false(trials_per_run, num_runs);
choice_lodiff_unc_onset_mask =  false(trials_per_run, num_runs);
choice_lodiff_con_onset_mask =  false(trials_per_run, num_runs);
choice_hidiff_unc_onset_mask =  false(trials_per_run, num_runs);
choice_hidiff_con_onset_mask =  false(trials_per_run, num_runs);
greencue_onset_mask = false(trials_per_run, num_runs);
cue_negfbktrial_onset_mask = false(trials_per_run, num_runs);
cue_lodiff_unc_onset_mask = false(trials_per_run, num_runs);
cue_lodiff_con_onset_mask = false(trials_per_run, num_runs);
cue_hidiff_unc_onset_mask = false(trials_per_run, num_runs);
cue_hidiff_con_onset_mask = false(trials_per_run, num_runs);
negfbk_lodiff_unc_onset_mask = false(trials_per_run, num_runs);
negfbk_lodiff_con_onset_mask = false(trials_per_run, num_runs);
negfbk_hidiff_unc_onset_mask = false(trials_per_run, num_runs);
negfbk_hidiff_con_onset_mask = false(trials_per_run, num_runs);
posfbk_lodiff_unc_onset_mask = false(trials_per_run, num_runs);
posfbk_lodiff_con_onset_mask = false(trials_per_run, num_runs);
posfbk_hidiff_unc_onset_mask = false(trials_per_run, num_runs);
posfbk_hidiff_con_onset_mask = false(trials_per_run, num_runs);
goal_posfbk_onset_mask = false(trials_per_run, num_runs);
goal_negfbk_onset_mask = false(trials_per_run, num_runs);
choice_miss_onset_mask = false(trials_per_run, num_runs);
choice_afternegfbk_miss_onset_mask = false(trials_per_run, num_runs);
greencue_miss_onset_mask = false(trials_per_run, num_runs);
cue_allneg_miss_onset_mask = false(trials_per_run, num_runs);
fbk_allneg_miss_onset_mask = false(trials_per_run, num_runs);


% open file (convert to ascii first)
if ~exist(fullfile(DATADIR,[DATAFILE_PRE subject DATAFILE_SUF]),'file')
    DATAFILE_PRE = 'PersTasknewcomp_fmriorder2_1234-';
    %error(sprintf('cannot find file %s',fullfile(DATADIR,[DATAFILE_PRE subject DATAFILE_SUF])));
end
if ~exist(fullfile(DATADIR,[DATAFILE_PRE subject DATAFILE_SUF]),'file')
    DATAFILE_PRE = 'PersTaskcumulative_fmri_1234-';
end
unicode2ascii(fullfile(DATADIR,[DATAFILE_PRE subject DATAFILE_SUF]),fullfile(DATADIR,[DATAFILE_PRE subject 'ascii' DATAFILE_SUF]));
datafid = fopen(fullfile(DATADIR,[DATAFILE_PRE subject 'ascii' DATAFILE_SUF]),'rt');

for irun = 1:num_runs
    trialnum = 1;
    last_offset = 4; %first screen at 4 seconds
    oneline = fgetl(datafid);
    while ischar(oneline)
        oneline = fgetl(datafid);
        if ~ischar(oneline) break; end
        if trialnum > trials_per_run
            break;
        end
        %get eventtype
        while (isempty(strfind(oneline,'EventType:'))) oneline = fgetl(datafid); end
        splitline = split(':',oneline); %1=uncontrol,2=control,3=poscue,4=choose,6=goalfbk
        event_type = str2num(strtok(splitline{2}));
        while (isempty(strfind(oneline,'Lose:'))) oneline = fgetl(datafid); end
        splitline = split(':',oneline);
        lose = str2num(strtok(splitline{2}));
        while (isempty(strfind(oneline,'ValueDiff:'))) oneline = fgetl(datafid); end
        splitline = split(':',oneline); %1=lodiff,2=hidiff
        valuediff = str2num(strtok(splitline{2}));
        while (isempty(strfind(oneline,'NewRound:'))) oneline = fgetl(datafid); end
        splitline = split(':',oneline);
        newround = str2num(strtok(splitline{2}));
        while (isempty(strfind(oneline,'PostCuejit:'))) oneline = fgetl(datafid); end
        splitline = split(':',oneline);
        postcuejit = str2num(strtok(splitline{2}));
        while (isempty(strfind(oneline,'PostFbkjit:'))) oneline = fgetl(datafid); end
        splitline = split(':',oneline);
        postfbkjit = str2num(strtok(splitline{2}));
        
        switch event_type
            case 4
                while (isempty(strfind(oneline,'ChoosePath.RT:'))) oneline = fgetl(datafid); end
                splitline = split(':',oneline);
                rt = str2num(strtok(splitline{2}));
                while (isempty(strfind(oneline,'postChoiceFixation.RT:'))) oneline = fgetl(datafid); end
                splitline = split(':',oneline);
                late_rt = str2num(strtok(splitline{2}));
                if (rt > 0)
                    rt_data(trialnum,irun) = rt;
                elseif (late_rt > 0)
                    rt_data(trialnum,irun) = 2000 + late_rt;
                else
                    rt_data(trialnum,irun) = NaN;
                end
                while (isempty(strfind(oneline,'Persist:'))) oneline = fgetl(datafid); end
                splitline = split(':',oneline);
                persist_choice = str2num(strtok(splitline{2}));
                if persist_choice == 1
                    persist_data(trialnum-1,irun) = 2;
                elseif persist_choice == 0
                    persist_data(trialnum-1,irun) = 1;
                end
                while (isempty(strfind(oneline,'Switch:'))) oneline = fgetl(datafid); end
                splitline = split(':',oneline);
                switch_choice = str2num(strtok(splitline{2}));
                if switch_choice == 1
                    switch_data(trialnum-1,irun) = 2;
                elseif switch_choice == 0
                    switch_data(trialnum-1,irun) = 1;
                end
                if isnan(rt_data(trialnum,irun))
                    choice_miss_onset_mask(trialnum, irun) = true;
                    if newround == 0 choice_afternegfbk_miss_onset_mask(trialnum, irun) = true; end
                elseif newround == 1
                    firstchoice_onset_mask(trialnum, irun) = true;
                end
                choice_onset_list(trialnum, irun) = last_offset;
                last_offset = last_offset + 2 + (postcuejit/1000);
                %choice onset masks get set based on the preceding
                %negfbk condition (cases 2 and 1 below)
            case 6
                while (isempty(strfind(oneline,'Goalimage:'))) oneline = fgetl(datafid); end
                splitline = split(':',oneline);
                if ~isempty(strfind(oneline,'yescap'))
                    goal_posfbk_onset_mask(trialnum,irun) = true;
                else
                    goal_negfbk_onset_mask(trialnum,irun) = true;
                end
                goal_onset_list(trialnum,irun) = last_offset;
                last_offset = last_offset + 4;
            case 3
                while (isempty(strfind(oneline,'greencue.RT:'))) oneline = fgetl(datafid); end
                splitline = split(':',oneline);
                rt = str2num(strtok(splitline{2}));
                cue_onset_list(trialnum,irun) = last_offset;
                if rt>0
                    greencue_onset_mask(trialnum,irun) = true;
                else
                    greencue_miss_onset_mask(trialnum,irun) = true;
                end
                last_offset = last_offset + 4 + (postfbkjit/1000);
            case 2
                while (isempty(strfind(oneline,'Cue.RT:'))) oneline = fgetl(datafid); end
                splitline = split(':',oneline);
                rt = str2num(strtok(splitline{2}));
                while (isempty(strfind(oneline,'postCueFixation.RT:'))) oneline = fgetl(datafid); end
                splitline = split(':',oneline);
                late_rt = str2num(strtok(splitline{2}));
                while (isempty(strfind(oneline,'*** LogFrame End ***'))) 
                    oneline = fgetl(datafid);
                    if ~isempty(strfind(oneline,'IncorrectKey: 1'))  %%count incorrectkey as a miss
                        rt = 0;
                        late_rt = 0;
                    end
                end
                if (rt > 0)
                    rt_data(trialnum,irun) = rt;
                elseif (late_rt > 0)
                    rt_data(trialnum,irun) = 2000 + late_rt;
                else
                    rt_data(trialnum,irun) = NaN;
                end
                cue_onset_list(trialnum,irun) = last_offset;
                fbk_onset_list(trialnum,irun) = last_offset + 2 + (postcuejit/1000);
                last_offset = last_offset + 2 + (postcuejit/1000) + 2 + (postfbkjit/1000);
                if isnan(rt_data(trialnum,irun))
                    cue_allneg_miss_onset_mask(trialnum,irun) = true;
                    fbk_allneg_miss_onset_mask(trialnum,irun) = true;
                    if lose == 1 cue_negfbktrial_onset_mask(trialnum,irun) = true; end
                elseif valuediff == 1 && lose == 1
                    cue_lodiff_con_onset_mask(trialnum,irun) = true;
                    negfbk_lodiff_con_onset_mask(trialnum,irun) = true;
                    cue_negfbktrial_onset_mask(trialnum,irun) = true;
                    choice_lodiff_con_onset_mask(trialnum+1,irun) = true; %set mask for following choice
                elseif valuediff == 2 && lose == 1
                    cue_hidiff_con_onset_mask(trialnum,irun) = true;
                    negfbk_hidiff_con_onset_mask(trialnum,irun) = true;
                    cue_negfbktrial_onset_mask(trialnum,irun) = true;
                    choice_hidiff_con_onset_mask(trialnum+1,irun) = true; %set mask for following choice
                elseif valuediff == 1 && lose == 0
                    cue_lodiff_con_onset_mask(trialnum,irun) = true;
                    posfbk_lodiff_con_onset_mask(trialnum,irun) = true;
                elseif valuediff == 2 && lose == 0
                    cue_hidiff_con_onset_mask(trialnum,irun) = true;
                    posfbk_hidiff_con_onset_mask(trialnum,irun) = true;
                end
            case 1
                while (isempty(strfind(oneline,'Cue.RT:'))) oneline = fgetl(datafid); end
                splitline = split(':',oneline);
                rt = str2num(strtok(splitline{2}));
                while (isempty(strfind(oneline,'postCueFixation.RT:'))) oneline = fgetl(datafid); end
                splitline = split(':',oneline);
                late_rt = str2num(strtok(splitline{2}));
                if (rt > 0)
                    rt_data(trialnum,irun) = rt;
                elseif (late_rt > 0)
                    rt_data(trialnum,irun) = 2000 + late_rt;
                else
                    rt_data(trialnum,irun) = NaN;
                end
                cue_onset_list(trialnum,irun) = last_offset;
                fbk_onset_list(trialnum,irun) = last_offset + 2 + (postcuejit/1000);
                last_offset = last_offset + 2 + (postcuejit/1000) + 2 + (postfbkjit/1000);
                if isnan(rt_data(trialnum,irun))
                    cue_allneg_miss_onset_mask(trialnum,irun) = true;
                    fbk_allneg_miss_onset_mask(trialnum,irun) = true;
                    if lose == 1 cue_negfbktrial_onset_mask(trialnum,irun) = true; end
                elseif valuediff == 1 && lose == 1
                    cue_lodiff_unc_onset_mask(trialnum,irun) = true;
                    negfbk_lodiff_unc_onset_mask(trialnum,irun) = true;
                    cue_negfbktrial_onset_mask(trialnum,irun) = true;
                    choice_lodiff_unc_onset_mask(trialnum+1,irun) = true; %set mask for following choice
                elseif valuediff == 2 && lose == 1
                    cue_hidiff_unc_onset_mask(trialnum,irun) = true;
                    negfbk_hidiff_unc_onset_mask(trialnum,irun) = true;
                    cue_negfbktrial_onset_mask(trialnum,irun) = true;
                    choice_hidiff_unc_onset_mask(trialnum+1,irun) = true; %set mask for following choice
                elseif valuediff == 1 && lose == 0
                    cue_lodiff_unc_onset_mask(trialnum,irun) = true;
                    posfbk_lodiff_unc_onset_mask(trialnum,irun) = true;
                elseif valuediff == 2 && lose == 0
                    cue_hidiff_unc_onset_mask(trialnum,irun) = true;
                    posfbk_hidiff_unc_onset_mask(trialnum,irun) = true;
                end
        end
        trialnum = trialnum + 1;
    end % while get line loop
end %for irun
fclose(datafid);

%now write the onset files, 1 per condition per run
% if ~exist(fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets'));
%     mkdir(fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets'));
% end

for irun = 1:num_runs
    if ~exist(fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun)));
        mkdir(fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun)));
    end
    %write onset in seconds, duration in seconds, weight=1 with tabs separating each
    %first choice of each round
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond001.txt');
    Mwrite = [ choice_onset_list(firstchoice_onset_mask(:,irun),irun) ...
        choice_duration_list(firstchoice_onset_mask(:,irun),irun) ...
        weights_list(firstchoice_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %choice after lodiff unc neg
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond002.txt');
    Mwrite = [ choice_onset_list(choice_lodiff_unc_onset_mask(:,irun)&~choice_miss_onset_mask(:,irun),irun) ...
        choice_duration_list(choice_lodiff_unc_onset_mask(:,irun)&~choice_miss_onset_mask(:,irun),irun) ...
        weights_list(choice_lodiff_unc_onset_mask(:,irun)&~choice_miss_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %choice after lodiff con neg
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond003.txt');
    Mwrite = [ choice_onset_list(choice_lodiff_con_onset_mask(:,irun)&~choice_miss_onset_mask(:,irun),irun) ...
        choice_duration_list(choice_lodiff_con_onset_mask(:,irun)&~choice_miss_onset_mask(:,irun),irun) ...
        weights_list(choice_lodiff_con_onset_mask(:,irun)&~choice_miss_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %choice after hidiff unc neg
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond004.txt');
    Mwrite = [ choice_onset_list(choice_hidiff_unc_onset_mask(:,irun)&~choice_miss_onset_mask(:,irun),irun) ...
        choice_duration_list(choice_hidiff_unc_onset_mask(:,irun)&~choice_miss_onset_mask(:,irun),irun) ...
        weights_list(choice_hidiff_unc_onset_mask(:,irun)&~choice_miss_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %choice after hidiff con neg
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond005.txt');
    Mwrite = [ choice_onset_list(choice_hidiff_con_onset_mask(:,irun)&~choice_miss_onset_mask(:,irun),irun) ...
        choice_duration_list(choice_hidiff_con_onset_mask(:,irun)&~choice_miss_onset_mask(:,irun),irun) ...
        weights_list(choice_hidiff_con_onset_mask(:,irun)&~choice_miss_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %positive cue and fbk
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond006.txt');
    Mwrite = [ cue_onset_list(greencue_onset_mask(:,irun),irun) ...
        cuefbk_duration_list(greencue_onset_mask(:,irun),irun) ...
        weights_list(greencue_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %cue lodiff unc
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond007.txt');
    Mwrite = [ cue_onset_list(cue_lodiff_unc_onset_mask(:,irun),irun) ...
        cue_duration_list(cue_lodiff_unc_onset_mask(:,irun),irun) ...
        weights_list(cue_lodiff_unc_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %cue lodiff con
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond008.txt');
    Mwrite = [ cue_onset_list(cue_lodiff_con_onset_mask(:,irun),irun) ...
        cue_duration_list(cue_lodiff_con_onset_mask(:,irun),irun) ...
        weights_list(cue_lodiff_con_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %cue hidiff unc
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond009.txt');
    Mwrite = [ cue_onset_list(cue_hidiff_unc_onset_mask(:,irun),irun) ...
        cue_duration_list(cue_hidiff_unc_onset_mask(:,irun),irun) ...
        weights_list(cue_hidiff_unc_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %cue hidiff con
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond010.txt');
    Mwrite = [ cue_onset_list(cue_hidiff_con_onset_mask(:,irun),irun) ...
        cue_duration_list(cue_hidiff_con_onset_mask(:,irun),irun) ...
        weights_list(cue_hidiff_con_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %negfbk lodiff unc
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond011.txt');
    Mwrite = [ fbk_onset_list(negfbk_lodiff_unc_onset_mask(:,irun),irun) ...
        fbk_duration_list(negfbk_lodiff_unc_onset_mask(:,irun),irun) ...
        weights_list(negfbk_lodiff_unc_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %negfbk lodiff con
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond012.txt');
    Mwrite = [ fbk_onset_list(negfbk_lodiff_con_onset_mask(:,irun),irun) ...
        fbk_duration_list(negfbk_lodiff_con_onset_mask(:,irun),irun) ...
        weights_list(negfbk_lodiff_con_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %negfbk hidiff unc
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond013.txt');
    Mwrite = [ fbk_onset_list(negfbk_hidiff_unc_onset_mask(:,irun),irun) ...
        fbk_duration_list(negfbk_hidiff_unc_onset_mask(:,irun),irun) ...
        weights_list(negfbk_hidiff_unc_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %negfbk hidiff con
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond014.txt');
    Mwrite = [ fbk_onset_list(negfbk_hidiff_con_onset_mask(:,irun),irun) ...
        fbk_duration_list(negfbk_hidiff_con_onset_mask(:,irun),irun) ...
        weights_list(negfbk_hidiff_con_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %posfbk lodiff unc
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond015.txt');
    Mwrite = [ fbk_onset_list(posfbk_lodiff_unc_onset_mask(:,irun),irun) ...
        fbk_duration_list(posfbk_lodiff_unc_onset_mask(:,irun),irun) ...
        weights_list(posfbk_lodiff_unc_onset_mask(:,irun),irun)];
    if isempty(Mwrite) Mwrite = [ 0 0 0]; end
    dlmwrite(outfilename,Mwrite,'\t');
    %posfbk lodiff con
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond016.txt');
    Mwrite = [ fbk_onset_list(posfbk_lodiff_con_onset_mask(:,irun),irun) ...
        fbk_duration_list(posfbk_lodiff_con_onset_mask(:,irun),irun) ...
        weights_list(posfbk_lodiff_con_onset_mask(:,irun),irun)];
    if isempty(Mwrite) Mwrite = [ 0 0 0]; end
    dlmwrite(outfilename,Mwrite,'\t');
    %posfbk hidiff unc
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond017.txt');
    Mwrite = [ fbk_onset_list(posfbk_hidiff_unc_onset_mask(:,irun),irun) ...
        fbk_duration_list(posfbk_hidiff_unc_onset_mask(:,irun),irun) ...
        weights_list(posfbk_hidiff_unc_onset_mask(:,irun),irun)];
    if isempty(Mwrite) Mwrite = [ 0 0 0]; end
    dlmwrite(outfilename,Mwrite,'\t');
    %posfbk hidiff con
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond018.txt');
    Mwrite = [ fbk_onset_list(posfbk_hidiff_con_onset_mask(:,irun),irun) ...
        fbk_duration_list(posfbk_hidiff_con_onset_mask(:,irun),irun) ...
        weights_list(posfbk_hidiff_con_onset_mask(:,irun),irun)];
    if isempty(Mwrite) Mwrite = [ 0 0 0]; end
    dlmwrite(outfilename,Mwrite,'\t');
    %goal pos
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond019.txt');
    Mwrite = [ goal_onset_list(goal_posfbk_onset_mask(:,irun),irun) ...
        goal_duration_list(goal_posfbk_onset_mask(:,irun),irun) ...
        weights_list(goal_posfbk_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %goal neg
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond020.txt');
    Mwrite = [ goal_onset_list(goal_negfbk_onset_mask(:,irun),irun) ...
        goal_duration_list(goal_negfbk_onset_mask(:,irun),irun) ...
        weights_list(goal_negfbk_onset_mask(:,irun),irun)];
    dlmwrite(outfilename,Mwrite,'\t');
    %misses
    %choice misses
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond021.txt');
    Mwrite = [ choice_onset_list(choice_miss_onset_mask(:, irun),irun) ...
        choice_duration_list(choice_miss_onset_mask(:, irun),irun) ...
        weights_list(choice_miss_onset_mask(:, irun),irun)];
    if isempty(Mwrite) Mwrite = [ 0 0 0]; end
    dlmwrite(outfilename,Mwrite,'\t');
    %fprintf('missed choices for subject %s, run %d: %d\n',subject,irun,length(weights_list(choice_miss_onset_mask(:, irun),irun)));
    fprintf('missed choices (after negfbk) for subject %s, run %d: %d\n',subject,irun,length(weights_list(choice_afternegfbk_miss_onset_mask(:, irun)&~firstchoice_onset_mask(:,irun),irun)));
    %poscue misses
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond022.txt');
    Mwrite = [ cue_onset_list(greencue_miss_onset_mask(:, irun),irun) ...
        cuefbk_duration_list(greencue_miss_onset_mask(:, irun),irun) ...
        weights_list(greencue_miss_onset_mask(:, irun),irun)];
    if isempty(Mwrite) Mwrite = [ 0 0 0]; end
    dlmwrite(outfilename,Mwrite,'\t');
    %fprintf('missed poscues for subject %s, run %d: %d\n',subject,irun,length(weights_list(greencue_miss_onset_mask(:, irun),irun)));
    %cue unc/con misses
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond023.txt');
    Mwrite = [ cue_onset_list(cue_allneg_miss_onset_mask(:, irun),irun) ...
        cue_duration_list(cue_allneg_miss_onset_mask(:, irun),irun) ...
        weights_list(cue_allneg_miss_onset_mask(:, irun),irun)];
    if isempty(Mwrite) Mwrite = [ 0 0 0]; end
    dlmwrite(outfilename,Mwrite,'\t');
    %fprintf('missed negcues (all) for subject %s, run %d: %d\n',subject,irun,length(weights_list(cue_allneg_miss_onset_mask(:, irun),irun)));
    %fprintf('missed negcues (fail trials) for subject %s, run %d: %d\n',subject,irun, ...
    %    length(weights_list(cue_allneg_miss_onset_mask(:, irun)&cue_negfbktrial_onset_mask(:,irun),irun)));
    %fbk unc/con misses
    outfilename = fullfile(EXP_DIR, sprintf('sub%s',subject), 'model', 'model001', 'onsets', sprintf('task001_run00%d',irun), 'cond024.txt');
    Mwrite = [ fbk_onset_list(fbk_allneg_miss_onset_mask(:, irun),irun) ...
        fbk_duration_list(fbk_allneg_miss_onset_mask(:, irun),irun) ...
        weights_list(fbk_allneg_miss_onset_mask(:, irun),irun)];
    if isempty(Mwrite) Mwrite = [ 0 0 0]; end
    dlmwrite(outfilename,Mwrite,'\t');
    %     if length(weights_list(irun,missed_event_mask(irun,:)))
    %         fprintf('missed resps for subject %s, run %d: %d\n',subject,irun,length(weights_list(irun,missed_event_mask(irun,:))));
    %     end
end % for irun = 1:nruns
% if WRITE_BEHAV
%     outfid = fopen (fullfile(DATADIR,'OptionLoss_behav.csv'),'at');
%     %write header in output file
% %     fprintf(outfid,'subject,');
% %     for irun = 1:files_to_read
% %         fprintf(outfid,'varstay_run%d,varswitch_run%d,consstay_run%d,consswitch_run%d,varchoice_noresp%d,conschoice_noresp%d,varshrinknoresp%d,consshrinknoresp%d,',irun,irun,irun,irun,irun,irun,irun,irun);
% %     end %for irun
% %     fprintf(outfid,'\n');
%     fprintf(outfid,'%s,',subject);
%     for irun = 1:files_to_read
%         %!!! the statement below doesn't work for more than one run -- fix it when
%         %necessary !!!
%         fprintf(outfid,sprintf('%d,%d,%d,%d,%d,%d,%d,%d,', ...
%             sum(stay_choice_mask & var_block_mask), ...
%             sum(switch_choice_mask & var_block_mask), ...
%             sum(stay_choice_mask & cons_block_mask), ...
%             sum(switch_choice_mask & cons_block_mask), ...
%             sum(noresp_list & var_block_mask & choicecue_mask), ...
%             sum(noresp_list & cons_block_mask & choicecue_mask), ...
%             sum(noresp_list & var_block_mask & shrinkcue_mask), ...
%             sum(noresp_list & cons_block_mask & shrinkcue_mask)));
%     end %for irun
%     fprintf(outfid,'\n');
%     fclose(outfid);
% end
% save(['subject' subject '_behavior']);

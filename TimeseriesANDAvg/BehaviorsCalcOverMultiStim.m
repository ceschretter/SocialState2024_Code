
%% This code is used for creating timeseries and averages following processing on the FlyDisco pipeline 
%Inputs - 
%Conditions: dirs containing expdirs for each experimental type. I
%sorted the expdirs by genotype.
%Behaviors: scores files from the classifier used
%Pre: number of frames that is the pre-stimulus period
%Post: number of frames that followed this period that you'd like to
%calculate over
%Num_stims: this is the number of stimuli you've provided
%Alignment_stim: this is the stimulus you'd like to be your reference
%Bin: number of frames you are binning by 
%FPS: frames per second for your camera
%Outputdir: this is where you'd like your excel file to go
%savefilestr: this is the name you'd like to give your excel file

%Output –
% This produces an excel file which has the percent of trajectories
% exhibiting the behavior per bin calculated over the timeframe you've
% given

%These conditions can be changed to the genotype names, Note: conditions cannot start with a number
conditions = {'Group1','Group2'};

%load in the behavior scores you'd like to look at
behaviors = {'scores_FeAggfRGBub.mat'}; %This is from the female aggression classifier used in the paper
behaviors_short = {'FeAgg'};

% set the following
writecsv = true;
plotfig = false;
pre = 2550;
post = 58650;
num_stims = 6;
alignment_stim = 1;
tot = pre+post;
post = post-1;
bin = 60;
fps = 169.99;
outputdir = '/groups/rubin/home/schretterc/Documents/FlyDiscoAnalysis_ExptsToAnalyze_Test/VisionProject_LCAnalysis/otherLCs';
savefilestr = ['PercentageofFlies_',behaviors{1}(1:end-4),'_','Meanbinsize',num2str(bin), '_',datestr(now,'yyyymmddTHHMMSS'),'.csv'];
savefilename = fullfile(outputdir,savefilestr);

%The tot needs to be evenly divisible by bin #
binlabels = nan(1,tot/bin);
binlabels(1) = mean(1:bin)*(1/169.99);
for i = 1:tot/bin-1
    binlabels(i+1) = mean(i*bin+1:i*bin+bin)*(1/fps);
end


if writecsv
    fid = fopen(savefilename,'w');

end
if plotfig
    figure;
end

for cond = 1:numel(conditions)
    condition = conditions{cond};

    switch condition
        %need to use the same naming as defined above
        case 'Group1'
            % Need to change this to each dir path would like to use
            datadirs = {'/Group1Dir'};
        case 'Group2'
            % Need to change this to each dir path would like to use
            datadirs = {'/Group2Dir'};
    end

    % make dir list
    datadir = datadirs{:};
    expdirs = dir(datadir);
    expdirs(ismember( {expdirs.name}, {'.', '..','.DS_Store'})) = [];

    combinedbowl_meanperframe = nan(numel(expdirs),tot);
    binned_combinedbowl_meanperframe  = nan(numel(expdirs),tot/bin);


    % write csv
    if writecsv
        fprintf(fid,'%s \n',condition);

    end
    indicatorstarts = nan(numel(expdirs),num_stims);
    for ee = 1:numel(expdirs)
        expdir = fullfile(datadir,expdirs(ee).name);
        % load indicator data
        indicatorfile = 'indicatordata.mat';
        if exist(fullfile(expdir,indicatorfile),'file')
            load(fullfile(expdir,indicatorfile),'indicatorLED')
        else
            sprintf('No indicatorLED file')
            return
        end
        %load in scores file
        scores = 'scores_FeAggfRGBub.mat';
        if exist(fullfile(expdir,scores),'file')
            load(fullfile(expdir,scores),'allScores')
        else
            sprintf('No trx file')
            return
        end

        nflies = numel(allScores.postprocessed);
        currbowl_data = nan(nflies,max(allScores.tEnd));
        for i = 1:nflies
            currbowl_data(i,allScores.tStart(i):allScores.tEnd(i)) = allScores.postprocessed{i}(allScores.tStart(i):allScores.tEnd(i));
        end
        % pull stimulus timing from all experiments
        indicatorstarts(ee,:) = indicatorLED.startframe;

        % find frame for data alignment based on alignment_stim
        stimstart = indicatorLED.startframe(alignment_stim);
        currbowl_data_cropped = currbowl_data(:,stimstart-pre:stimstart+post);
        for i = 1:size(currbowl_data_cropped,2)
            currbowl_percbehavior(i) = (nansum((currbowl_data_cropped(:,i)))/numel(find(~isnan(currbowl_data_cropped(:,i)))).*100);
        end
        combinedbowl_percbehavior(ee,:) = currbowl_percbehavior;
        binned_combinedbowl_percbehavior(ee,1) = mean(currbowl_percbehavior(1:bin));
        binlabels(1) = mean(1:bin)*(1/169.99);
        for i = 1:tot/bin-1
            binned_combinedbowl_percbehavior(ee,i+1) = mean(currbowl_percbehavior(i*bin+1:i*bin+bin));

        end

        % write data to csv
        if writecsv

            for i = 1:numel(binned_combinedbowl_percbehavior(ee,:))
                fprintf(fid, '%.3f, ',binned_combinedbowl_percbehavior(ee,i));
            end
            fprintf(fid, '\n');
        end
    end

    % test that indicators are evenly spaced.
    diff(indicatorLED.startframe)

    % check indicator timing matches – if your indicator data is off then
    % this will fail
    a = num2cell(diff(indicatorstarts'));
    assert(isequal(a{:}),'indicator timing is not equal')


    %Before progressing you can plot your data if you selected plotfig = true
    if plotfig
        subplot(2,2,cond)
        patch([pre/bin,(pre+900)/bin,(pre+900)/bin,pre/bin],[0, 0, 30 , 30],'r','facealpha',.25)
        hold on
        plot(mean(binned_combinedbowl_percbehavior,1),'r')
        title(condition,'Interpreter','none')
    end

    % create binned inidcator data and adjust timestamps using last experiment loaded for each
    % condition
    indicatordigital_cropped = indicatorLED.indicatordigital(stimstart-pre:stimstart+post);

    binned_indicatordata = nan(1,size(binned_combinedbowl_percbehavior,2));
    binned_indicatordata(1) = mean(indicatordigital_cropped(1:bin));
    for i = 1:tot/bin-1
        binned_indicatordata(i+1) = mean(indicatordigital_cropped(i*bin+1:i*bin+bin));
    end
    binned_indicatordata = binned_indicatordata >= .5; %set threshold to greater than 50% of the values are a 1

    % align timestamps to movie
    offset = (stimstart-pre)*(1/fps);
    binlabels = binlabels+offset;

    if writecsv
        % add binned indicator data
        fprintf(fid,'binned digital indicator \n');
        for i = 1:tot/bin
            fprintf(fid,'%f, ',binned_indicatordata(i));
        end
        fprintf(fid,'\n');

        % add timestamps aligned to stimuli
        fprintf(fid,'timestamps \n');
        for i = 1:tot/bin
            fprintf(fid,'%.2f, ',binlabels(i));
        end
        fprintf(fid,'\n');

    end
end

%finished writing excel file
if writecsv
    fclose(fid);
end


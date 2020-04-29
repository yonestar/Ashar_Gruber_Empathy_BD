clear all; close all
basedir = '~/Repositories/Ashar_Gruber_Empathy_BD';

load(fullfile(basedir, 'project_workspace.mat'))

% patients are group 1; healthies are group 2

%% IF LOADED IN DATA IN ABOVE CELL, AND DON'T NEED TO MAKE CHANGES TO DATA, 
%  CAN SKIP TO LINE 274

%%

newstruct = struct('happy1',[],'sad1',[],'tender1',[], 'happy2',[],'sad2',[],'tender2',[]);
%[amy, angela, anne, anthony, bill, charles, crystal, destiny, fanny, isaiah, jasmine, jesse, jessica, jodi, kelly, lance, mike, natalie, nathan, robert, ron, valarie, wallace, will] = deal(newstruct);

%% load task data

tabl = readtable(fullfile(basedir,'BD_empathic_emotions_continuous_ratings__Survey_C_Part_1_Final 2.csv'));
tabl(1,:) = [];  % first row is dummy row

%% show all duplicates in task data.  DONT NEED TO RE_RUN
if 0
    for i=1:height(tabl)

        tmp = strcmp(tabl.ID(i), tabl.ID);
        if sum(tmp) > 1
            fprintf('Duplicate for ID %s\n', tabl.ID{i});
            tabl(tmp,:)
        end
    end
end
%% DROP DUPLICATES. RUN EVERY TIME.  
%for 135_Y, the first row is the good one
tmp = find(strcmp('135_Y', tabl.ID))
to_drop = tmp(2:3);
to_drop = to_drop';

% 175_Y, 171_B 245_Y 248_Y 263_Y and more, the first one is good
tmp = find(strcmp('175_Y', tabl.ID))
to_drop = [to_drop tmp(2)];

tmp = find(strcmp('171_B', tabl.ID))
to_drop = [to_drop tmp(2)];

tmp = find(strcmp('245_Y', tabl.ID))
to_drop = [to_drop tmp(2)];

tmp = find(strcmp('248_Y', tabl.ID))
to_drop = [to_drop tmp(2)];

tmp = find(strcmp('263_Y', tabl.ID))
to_drop = [to_drop tmp(2)];

tmp = find(strcmp('69_Y', tabl.ID))
to_drop = [to_drop tmp(2)];

%tmp = find(strcmp('173_Y', tabl.ID))
%to_drop = [to_drop tmp(2)];

tmp = find(strcmp('183_Y', tabl.ID))
to_drop = [to_drop tmp(2)];

tmp = find(strcmp('374_Y', tabl.ID))
to_drop = [to_drop tmp(2)];

% 55_Y 321_Y 347_Y 49_Y  2nd is good
tmp = find(strcmp('55_Y', tabl.ID))
to_drop = [to_drop tmp(1)];

tmp = find(strcmp('321_Y', tabl.ID))
to_drop = [to_drop tmp(1)];

tmp = find(strcmp('347_Y', tabl.ID))
to_drop = [to_drop tmp(1)];

tmp = find(strcmp('49_Y', tabl.ID))
to_drop = [to_drop tmp(1)];

% 65_Y both bad
tmp = find(strcmp('65_Y', tabl.ID))
to_drop = [to_drop tmp'];

tabl(to_drop,:) = [];

%% check I got them all. yes.
length(unique(tabl.ID)) == height(tabl)

%% remove Ps w/ no task data 

% ppl missing data for 1st and 2nd movie.  
nodata = (cellfun(@isempty,tabl.movie1_valence_vector)) & (cellfun(@isempty,tabl.movie2_valence_vector)) & (cellfun(@isempty,tabl.movie3_valence_vector)) & (cellfun(@isempty,tabl.movie4_valence_vector));

fprintf('Removing %d participants who have no task data (did not do the study)\n', sum(nodata));
tabl(nodata,:) = [];

%% remove one ineligible P, ID = 201 (see email from June 4/23/2020)

todrop = find(~cellfun(@isempty, strfind(tabl.ID, '201')));
tabl(todrop,:) = [];

%% Load survey data
surv = readtable(fullfile(basedir, 'Additional+Measures+for+Empathic+Emotion+Continuous+Rating+Task+-++Survey+C%2C+Part+2+Final_February+19%2C+2018_10.39.csv'));

% pull in ASRM and BDI-SF -- symptom severity
wh_ss = find(contains(surv.Properties.VariableNames, 'ASRM'));
wh_ss2 = find(contains(surv.Properties.VariableNames, 'BDI'));
wh_ss3 = find(contains(surv.Properties.VariableNames, 'IRI'));

surv.Properties.VariableNames(wh_ss)
surv.Properties.VariableNames(wh_ss2)
surv.Properties.VariableNames(wh_ss3)

% compute total scores
surv.ASRMtot = mean(surv{:, wh_ss},2);
surv.BDISFtot = mean(surv{:, wh_ss2},2);
surv.IRItot = mean(surv{:, wh_ss3}, 2);


mysurv = surv(:, {'ID', 'ASRMtot', 'BDISFtot', 'IRItot', 'SES', 'Income', 'Employment', 'Education', 'Gender', 'Age', 'Ethnicity' });

%% merge surveys w/ other data

% no one has a survey but no task data
numel(setdiff(mysurv.ID, tabl.ID))

% 10 ppl have task data but not survey
numel(setdiff(tabl.ID, mysurv.ID)) % returns the values in A that are not in B

tabl = outerjoin(tabl, mysurv, 'Keys', 'ID', 'MergeKeys', 1); % this will add NaN for the Ps without survey data

%% load group assignments
group = readtable(fullfile(basedir, 'group_assignments.csv'));
group.Properties.VariableNames{1} = 'ID';

%% in group, but no task data -- no one :)
setdiff(group.ID, tabl.ID)

% have task data, but no group data. no one :)
setdiff(tabl.ID, group.ID)

%% merge
tabl = join(tabl, group, 'Keys', 'ID');

%% how do symtpoms related to group status?
create_figure('s by g', 2,3);
boxplot(tabl.ASRMtot, tabl.group); title('ASRM'), xlabel('Group')
subplot(2,3,2)
boxplot(tabl.BDISFtot, tabl.group); title('BDI-SF'), xlabel('Group')

[~,p, ~, stats] = ttest2(tabl.ASRMtot(tabl.group==1),tabl.ASRMtot(tabl.group==2))
[~,p, ~, stats] = ttest2(tabl.BDISFtot(tabl.group==1),tabl.BDISFtot(tabl.group==2))

subplot(2,3,4); hist(tabl.ASRMtot); title('ASRM');
subplot(2,3,5); hist(tabl.BDISFtot); title('BDI-SF');

% corr of BDI and ASRM
subplot(2,3,3); scatter(tabl.BDISFtot(tabl.group==1), tabl.ASRMtot(tabl.group==1)); lsline; title('Group 1'); xlabel('BDI-SF'); ylabel('ASRM')
subplot(2,3,6); scatter(tabl.BDISFtot(tabl.group==2), tabl.ASRMtot(tabl.group==2)); lsline; title('Group 2'); xlabel('BDI-SF'); ylabel('ASRM')

%r = nancorr(tabl.BDISFtot, tabl.ASRMtot)


%% remove Ps who failed Careless Responder Q's

% check this subject's data quality
% careless respondng questions
wh_drop = ~strcmp(tabl.Q29_1,'') & ~strcmp(tabl.Q29_1,'2');
fprintf('Dropping %d Ss who failed careless responding question 29_1\n', sum(wh_drop));

wh_drop = ~strcmp(tabl.Q77,'') & ~strcmp(tabl.Q77,'4');
fprintf('Dropping %d Ss who failed careless responding question 77\n', sum(wh_drop));

tabl = tabl(~wh_drop, :);

%% remove Ps who's task data is all 0 -- presume error, not real data

wh_drop = zeros(height(tabl),1);
for i=1:height(tabl)     

    % check that didn't say all 0's for every valence_vector
    ok = 0;
    for m = 1:24
        vvname = ['movie' num2str(m) '_valence_vector'];
        vv = tabl{i, vvname}{1};
        vv = str2num(strrep(vv, '"', ''));
        if any(vv)
            ok = 1;
        end
    end
    
    if ~ok
        wh_drop(i) = 1;
    end
end

fprintf('Dropping %d Ss who had all zeros in every valence vector\n', sum(wh_drop));
        
tabl = tabl(~wh_drop, :);

%% Drop Ss with almost no data or implausible data.  See code below demonstrating this.  
% dropping them up here, early in code, so their data gets dropped from all
% data structures

wh_drop = [23 44 55 28];
tabl(wh_drop,:) = [];


%% --- ALL DATA IN TABL IS CONSIDERED 'GOOD CLEAN DATA' FROM HERE ON OUT ---

%% preproc and add data to structs
qualtrics = struct();
subjdat = {};

%create_figure('tvs', 3, 4)
for i=1:height(tabl)     
    
    % save data into preproc'ed format
    tvs{i} = [];
    for m = 1:24 % 24 movies
    
        % set up var names
        mname = ['movie_name' num2str(m)];
        emoname = ['emotion' num2str(m)];
        vvname = ['movie' num2str(m) '_valence_vector'];
        tvname = ['movie' num2str(m) '_time_vector'];
        
        % extract data for this row for this movie
        targ = tabl{i, mname}{1};
        emo = tabl{i, emoname}{1};
        vv = tabl{i, vvname}{1};
        tv = tabl{i, tvname}{1};
        
        % add group to emo
        emo = [emo num2str(tabl.group(i))];
        
        % if no data provided for this trial (i.e. dropped out early), skip
        if isempty(vv), continue, end
        
        % remove quotatation marks and cast to double; remove spaces in
        % emonames
        vv = str2num(strrep(vv, '"', ''));
        tv = str2num(strrep(tv, '"', ''));
        emo = strrep(emo, ' ', '_');
        
        tvs{i} = [tvs{i}; diff(tv)'];  %aggregate all delays between ratings
        
        % make sure tv (time vector) does not have repeat values.  this
        % should be quite rare, presumably caused by a rare software
        % glitch.  if repeat values are present, drop them
        [~,tmp]=unique(tv); tmp = setdiff(1:length(tv),tmp); 
        if ~isempty(tmp)
            tv(tmp) = []; vv(tmp) = [];
            fprintf('Dropped a repeated timepoint for subj %d movie %d.  Not a concern.\n', i, m);
        end
        
        %interp1 valence vector based on time vector to every .5 seconds
        vv_interp = interp1(tv, vv, 0:500:33000) ;
        
        % see if interpolation looks right
        %create_figure('testinterp'); plot(tv, vv, 'b'); plot(0:500:33000, vv_interp, 'g');
        
        % save in QUALTRICS data structure, grouped by targ
        %append interpolated valence vector to targ.emo
        if ~isfield(qualtrics, targ), qualtrics.(targ) = newstruct(); end
        qualtrics.(targ).(emo) = [qualtrics.(targ).(emo); vv_interp]; 
        
        % save in SUBJDAT data structure, grouped by subj
        if length(subjdat)<i % first trial for this subj
            subjdat{i} = struct('happy',[],'sad',[],'tender',[]);
        end
        % append it
        subjdat{i}.(emo(1:end-1)) = [ subjdat{i}.(emo(1:end-1)); vv_interp];
        
    end
    
    if sum(tvs{i} > 1100) > 10
        fprintf('Row %d had %d delays between recordings > 1100ms across all trials\n', i, sum(tvs{i} > 1100));
    end
    %subplot(3,4,i-60); hist(tvs)
end


%% SAVE RESULTS SO DON'T NEED TO REPEAT THE ABOVE
save(fullfile(basedir, 'project_workspace.mat'))


%% visualize results, target by target

targs = fieldnames(qualtrics);% {'amy','angela','anne','anthony','bill','charles','crystal','destiny','fanny','isaiah','jasmine','jesse','jessica','jodi','kelly','lance','mike','natalie','nathan','robert','ron','valarie','wallace','will'};
emos = fieldnames(qualtrics.angela);


create_figure('12 targs', 3, 4); 

% set up colors. modified from pexp_plots.m to match up, b/c emos are in
% different order here
emos = {'happy1','sad1','tender1', 'happy2','sad2','tender2'};
mycolors = {[.9 .9 0] [0 0 1] [1 0 1] [.9 .9 0] [0 0 1] [1 0 1]};
linstyle = {'-' '-' '-' '-.' '-.' '-.'};
    
for i = 1:12
    thistarg = [];
    
    % to see second 12 targs
    t = i; % + 12;
    
    for e=1:length(emos)
        thistarg = [thistarg; nanmean(qualtrics.(targs{t}).(emos{e}))];
    end
    
    thistarg=resample(thistarg', 26, 67);
    
    subplot(3,4,i)
    for j=1:length(emos)
        plot(thistarg(:,j), 'LineWidth', 4, 'color', mycolors{j}, 'LineStyle', linstyle{j})
    end
    xlim([0 27])
    ylim([- 5 90])
    title(targs{t})
end
legend(emos)


%% visualize results, averaged across all targets, by group

create_figure('all targs by group'); 
mycolors = {[.9 .9 0] [0 0 1] [1 0 1] [.9 .9 0] [0 0 1] [1 0 1]};
linstyle = {'-' '-' '-' ':' ':' ':'};    

targs = fieldnames(qualtrics);
emos = {'happy1','sad1','tender1', 'happy2','sad2','tender2'};
    
totalemo = [];
for e=1:length(emos)
    thisemo = [];
        
    for i = 1:24    
        t = i;
        thisemo = [thisemo; qualtrics.(targs{t}).(emos{e})];
    end
    
    totalemo.(emos{e}) = thisemo';
    
    h(e) = lineplot_columns(thisemo, 'w', 3, 'color', mycolors{e}, 'shade', 'linestyle', linstyle{e}, 'CIs', 'marker', 'none');
    h2(e) = h(e).line_han;
end

%% 
legend(h2, {'happy', 'sad', 'tender'}, 'Location', 'best')
set(gca, 'XTickLabel', 0:5:35)
xlabel('Time (s)')
ylabel('Emotion intensity')
title('Group average emotional responses')

saveas(gcf, fullfile(basedir, '..', 'Writing', 'figures', 'Fig1.jpg'))

%% visual results per subject
mycolors = {[.9 .9 0] [0 0 1] [1 0 1] [.9 .9 0] [0 0 1] [1 0 1]};
myemos = {'happy' 'sad' 'tender'};

create_figure('Group 1, by subj', 6, 10);

for i=1:length(subjdat)

    subplot(6, 10, i);
    hold on
    
    for j=1:3
        dat = subjdat{i}.(myemos{j});
        sz(i,j) = size(dat, 1);
        
        if sz(i,j)==0
            continue; 
        else
            %lineplot_columns(dat, 'w', 3, 'color', mycolors{j}, 'shade', 'CIs', 'marker', 'none');
            plot(dat',  'color', mycolors{j})
        end
        
    end
  
    titl = sprintf('#%d, n=%d %d %d', i, sz(i,:));
    title(titl, 'FontSize', 12)
    ylim([-10 110]), xlim([0 70]);
end

% Conclusion:  drop #44 and #55, have almost no data.  Drop #23, data does
% not look real -- all flat lines at different intercepts.  Now also drop
% #28, which is the only subj who doesn't have at least 1 timecourse for
% all 3 emos.  I DO THIS DROP ABOVE, SO THIS DATA ALSO GETS DROPPED FROM
% THE QUALTRICS STRUCT.  


%% ml-glm fitting a 2nd order polynomial for each first level
X=[];
X(:,1) = [1:67]';
X(:,2) = X(:,1) .^ 2;
%X(:,3) = X(:,1) .^ 3;
%X(:,4) = X(:,1) .^ 4;

X = repmat({X}, size(subjdat));

% choose which emo. find each P's mean timecourse for that emo across bios.
emo='happy'%'tender';
Y = cellfun( (@(s) s.(emo)), subjdat, 'UniformOutput',0);
Y = cellfun( @nanmean, Y, 'UniformOutput',0);
Y = cellfun( @transpose, Y, 'UniformOutput',0);

% group parameter for 2nd level
seclev_grop = tabl.group;
seclev_grop(seclev_grop==1)=-1;
seclev_grop(seclev_grop==2)=1;

fprintf('------------- %s -------------\n', emo);
glmfit_multilevel(Y, X, seclev_grop, 'noplots', 'names', {'Intcpt', 'linear', 'quad'});

% to look at BDI-SF as moderator
glmfit_multilevel(Y, X, tabl.BDISFtot, 'noplots', 'names', {'Intcpt', 'linear', 'quad'});




%% load donation data

dons_tabl = readtable(fullfile(basedir, 'donation_items_only_cleaned.csv'));
figure; imagesc(table2array(dons_tabl(:, 2:end))); colorbar

dons_tabl.avgdon = nanmean(table2array(dons_tabl(:, 2:end)), 2);
tabl = join(tabl, dons_tabl, 'Keys', 'ID');

tabl.avgdon % no NaNs here -- OK.

tabl.avgdon_cents = (tabl.avgdon-1) * 10;


%% average donation per group
g=2;
min(tabl.avgdon_cents), max(tabl.avgdon_cents)
mean(tabl.avgdon_cents(tabl.group==g))
std(tabl.avgdon_cents(tabl.group==g))


%% group differences in donations? no

boxplot(tabl.avgdon, tabl.group); title('Average donation'), xlabel('Group'), 
set(gca, 'FontSize', 24, 'YTick', 1:2:11, 'YTickLabel', {'$0', '$0.20', '$0.40', '$0.60', '$0.80', '$1.00'})
[~,p, ~, stats] = ttest2(tabl.avgdon(tabl.group==1),tabl.avgdon(tabl.group==2))


%% group differences in SES? 

create_figure('SES', 1, 3);
boxplot(tabl.SES, tabl.group); title('SES Ladder'), xlabel('Group');
subplot(1,3,2);
boxplot(tabl.Income, tabl.group); title('Income'), xlabel('Group');
subplot(1,3,3);
boxplot(tabl.Education, tabl.group); title('Education'), xlabel('Group');
%set(gca, 'FontSize', 24);

[~,p] = ttest2(tabl.SES(tabl.group==1),tabl.SES(tabl.group==2))
[~,p] = ttest2(tabl.Income(tabl.group==1),tabl.Income(tabl.group==2))
[~,p] = ttest2(tabl.Education(tabl.group==1),tabl.Education(tabl.group==2))

% employment status 
histc(tabl.Employment(tabl.group==1), 1:5)'
histc(tabl.Employment(tabl.group==2), 1:5)'


%% for demographics

demographics = tabl(:, {'ID' 'Gender' 'SES' 'group' 'Ethnicity' 'Age' 'Employment' 'Education' 'Income' 'ASRMtot' 'BDISFtot' 'IRItot'})

writetable(demographics, fullfile(basedir, 'demographics_final_sample.csv'))
%% Run script for tutorial examples

%%
cfg = [];
cfg.layout = 'neuromag306mag.lay';
cfg.showcomment = 'no';
ft_multiplotER(cfg, d1{:}); title('Condition 1')

ft_multiplotER(cfg, d2{:}); title('Condition 2')

%%
cfg = [];
cfg.layout = 'elec1010.lay';
ft_multiplotER(cfg, d1{:})

%% 
ga1 = ft_timelockgrandaverage([], d1{:});
ga2 = ft_timelockgrandaverage([], d2{:});


%%
cfg = [];
cfg.layout = 'neuromag306mag.lay';
cfg.showcomment = 'no';
ft_multiplotER(cfg, ga1, ga2)

%%
cfg = [];
cfg.layout = 'elec1010.lay';
cfg.showcomment = 'no';
ft_multiplotER(cfg, ga1, ga2)

%% Combine 
ga0 = ft_timelockgrandaverage([], d1{:}, d2{:});

cfg = [];
cfg.layout = 'neuromag306mag.lay';
cfg.showcomment = 'no';
ft_multiplotER(cfg, ga0)

%% Find peak channel
t_win = [0.140 0.180];

cfg = [];
cfg.channel = 'megmag';
ga0 = ft_selectdata(cfg, ga0);

% Select windows
toi = ga0.time > t_win(1) & ga0.time < t_win(2);

% Find peak channel
[mx1, idx1] = max(mean(ga0.avg(:,toi),2));
pk_chan = ga0.label(idx1);


%% Select data
cfg = [];
cfg.channel     = pk_chan;
cfg.latency     = t_win;
cfg.avgovertime = 'yes';

tmp_d1 = cell(size(d1));
tmp_d2 = cell(size(d2));

for ii = 1:20
    tmp_d1{ii} = ft_selectdata(cfg, d1{ii});
    tmp_d2{ii} = ft_selectdata(cfg, d2{ii});
end


%% Peak
tmp_d1 = [tmp_d1{:}];
tmp_d2 = [tmp_d2{:}];

pk_d1 = [tmp_d1.avg];
pk_d2 = [tmp_d2.avg];

%% Plot
figure;
scatter(ones(size(pk_d1)).*randi([85,115], 1,length(pk_d1))/100, pk_d1, 'b'); hold on
scatter(ones(size(pk_d2)).*randi([185,215], 1,length(pk_d1))/100, pk_d2, 'r')
line([0.75, 1.25], [mean(pk_d1), mean(pk_d1)], 'Color','b')
line([1.75, 2.25], [mean(pk_d2), mean(pk_d2)], 'Color','r')
xlim([0.5 2.5])

%% 
[mean(pk_d1), std(pk_d1)]
[mean(pk_d2), std(pk_d2)]

[H, P, CI, t] = ttest(pk_d1, pk_d2)

%% Select single channel data
cfg = [];
cfg.channel     = pk_chan;
cfg.avgovertime = 'no';

chan_d1 = cell(size(d1));
chan_d2 = cell(size(d2));

for ii = 1:20
    chan_d1{ii} = ft_selectdata(cfg, d1{ii});
    chan_d2{ii} = ft_selectdata(cfg, d2{ii});
end

%% Caluclate grand average
ga_d1 = ft_timelockgrandaverage([], chan_d1{:});
ga_d2 = ft_timelockgrandaverage([], chan_d2{:});

%% Plot grqand averages
cfg = [];
cfg.channel = pk_chan;
cfg.showcomment = 'no';
ft_singleplotER(cfg, ga_d1, ga_d2)

%% Cluster-based permutation tests on single channel data
cfg = [];
cfg.method              = 'montecarlo';                 % Permutation
cfg.correctm            = 'cluster';                    % Cluster-based correction
cfg.statistic           = 'ft_statfun_depsamplesT';     % Dependent samples (i.e. paired data)
cfg.alpha               = 0.025;                        % Two tailes (0.05/2) 
cfg.tail                = 0;
cfg.clusterstatistic    = 'maxsum';
cfg.clusteralpha        = 0.05;
cfg.clustertail         = 0;
cfg.numrandomization    = 1000;

% Desing matrix
tmp_ivar  = [ones(1, length(d1)), ones(1, length(d2))*2];   % Indepentend variable (condition)
tmp_unit = [1:size(d1), 1:size(d2)];                        % Unit of observation (subject)
design = [tmp_ivar; tmp_unit]';                             % Transpose to make it "long" format

% imagesc(design)

cfg.design              = design;
cfg.ivar                = 1;
cfg.uvar                = 2;

stat_single = ft_timelockstatistics(cfg, chan_d1{:}, chan_d2{:});

%% Plot uncorrected
figure; hold on
plot(stat_single.time, stat_single.stat)
ylabel('t statistic'); xlabel('time')
line([stat_single.time(1) stat_single.time(end)], [tinv(0.05/500, 19), tinv(0.05/500, 19)], 'Color','b')
line([stat_single.time(1) stat_single.time(end)], [tinv(0.95, 19), tinv(0.95, 19)], 'Color','r')

% Bonferroni
line([stat_single.time(1) stat_single.time(end)], [tinv(0.025/500, 19), tinv(0.025/500, 19)], 'Color','b')
line([stat_single.time(1) stat_single.time(end)], [tinv(1-0.025/500, 19), tinv(1-0.025/500, 19)], 'Color','b')

%% Plot stat
% Copy stat mask
ga_d1.mask = stat_single.mask;

% Plot
cfg = [];
cfg.channel         = pk_chan;
cfg.maskparameter   = 'mask';
cfg.maskfacealpha   = 0.5;
cfg.showcomment = 'no'
ft_singleplotER(cfg, ga_d1, ga_d2)

figure;
subplot(2,1,1); 
hist(stats.negdistribution, 50)
crit = prctile(stats.negdistribution, 2.5);
line([crit, crit], [0 300], 'LineWidth', 2, 'Color', 'r');
line([stats.negclusters(1).clusterstat, stats.negclusters(1).clusterstat], [0 300], 'LineWidth', 2, 'Color', 'b', 'LineStyle','--');
title('Negative cluster T')

subplot(2,1,2); hist(stats.posdistribution, 100)
crit = prctile(stats.posdistribution, 97.5);
line([crit, crit], [0 300], 'LineWidth', 2, 'Color', 'r');
line([stats.posclusters(1).clusterstat, stats.posclusters(1).clusterstat], [0 300], 'LineWidth', 2, 'Color', 'b', 'LineStyle','--');
title('Positive cluster T')


%% Prepare neighbours
cfg = [];
cfg.method  = 'template';
cfg.layout  = 'neuromag306mag.lay';
cfg.channel = 'megmag';
neighbours_mag = ft_prepare_neighbours(cfg, d1{1});

cfg = [];
cfg.method  = 'template';
cfg.layout  = 'elec1010.lay';
cfg.channel = 'eeg';
neighbours_eeg = ft_prepare_neighbours(cfg, d1{1});

%% Plot neighbours for inspection
cfg = [];
cfg.neighbours = neighbours_mag;
cfg.senstype = 'meg';
ft_neighbourplot(cfg, d1{1});

%% 
cfg = [];
cfg.operation = 'x1-x2';
cfg.parameter = 'avg';
dif = cell(size(d1));
for ii = 1:nn
    dif{ii} = ft_math(cfg, d1{ii}, d2{ii});
end

dif_ga = ft_timelockgrandaverage([], dif{:});

%%
cfg = [];
cfg.layout = 'neuromag306mag.lay';
cfg.showcomment = 'no';
ft_multiplotER(cfg, dif{:})
ft_multiplotER(cfg, dif_ga)

cfg = [];
cfg.layout = 'elec1010.lay';
cfg.showlabels = 'yes';
cfg.showcomment = 'no';
ft_multiplotER(cfg, dif{:})

%% 
cfg = [];
cfg.method              = 'montecarlo'; 
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.clusteralpha        = 0.05;        
cfg.clusterstatistic    = 'maxsum';
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.025;
cfg.channel             = 'megmag';
cfg.numrandomization    = 1000;
cfg.neighbours          = neighbours_mag;

% Desing matrix
tmp_ivar  = [ones(1, length(d1)), ones(1, length(d2))*2];
tmp_unit = [1:size(d1), 1:size(d2)];
design = [tmp_ivar; tmp_unit]';             % Transpose to make it "long" format

cfg.design = design;
cfg.ivar = 1;
cfg.uvar = 2; 

stats = ft_timelockstatistics(cfg, d1{:}, d2{:});

imagesc(design)

%%
figure;

subplot(2,1,1); 
hist(stats.negdistribution, 50)
crit = prctile(stats.negdistribution, 2.5);
line([crit, crit], [0 80], 'LineWidth', 2, 'Color', 'r');
line([stats.negclusters(1).clusterstat, stats.negclusters(1).clusterstat], [0 80], 'LineWidth', 2, 'Color', 'b', 'LineStyle','--');
title('Negative cluster T')

subplot(2,1,2); hist(stats.posdistribution, 100)
crit = prctile(stats.posdistribution, 97.5);
line([crit, crit], [0 80], 'LineWidth', 2, 'Color', 'r');
line([stats.posclusters(1).clusterstat, stats.posclusters(1).clusterstat], [0 80], 'LineWidth', 2, 'Color', 'b', 'LineStyle','--');
title('Positive cluster T')

%% Plot
% copy stat mask
dif_ga.mask = stats.mask;

% Plot
figure;
cfg = [];
cfg.channel         = 'megmag';
cfg.maskparameter   = 'mask';
cfg.maskfacealpha   = 0.5;
cfg.layout          = 'neuromag306mag';
ft_multiplotER(cfg, dif_ga)

% copy stat mask
ga1.mask = stats.mask;

% Plot
figure;
cfg = [];
cfg.channel         = 'megmag';
cfg.maskparameter   = 'mask';
cfg.maskfacealpha   = 0.5;
cfg.layout          = 'neuromag306mag';
cfg.showcomment = 'no';
ft_multiplotER(cfg, ga1, ga2)

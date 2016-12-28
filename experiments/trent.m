clear
addpath('~/src/TRENTOOL3')
addpath('~/src/fieldtrip-20160522');
ft_defaults;

fs = 6103.52; % Hz
segs = loadcached('data/ARN038_segs2','segs');
segs = cellfun(@(s) s', segs, 'Un',0);
tt = 0:1/fs:(size(segs{1},2)-1)/fs;

data.fsample = fs;
data.trial = segs;
[data.time{1:numel(segs)}] = deal(tt);
data.label = {'2' '4' '6' '8'  '10' '12' '14' '16'};

cfgTEP = [];
cfgTEP.toi = tt([1 end]);
cfgTEP.channel = data.label;

% optimizing embedding
cfgTEP.optimizemethod ='ragwitz';  % criterion used
cfgTEP.ragdim         = 2:9;       % criterion dimension
cfgTEP.ragtaurange    = [0.2 0.4]; % range for tau
cfgTEP.ragtausteps    = 5;        % steps for ragwitz tau steps
cfgTEP.repPred        = fix(numel(tt)*3/4);  % size(data.trial{1,1},2)*(3/4); (was 100)

% scanning of interaction delays u
cfgTEP.predicttime_u = 20;      % minimum u to be scanned

% estimator
cfgTEP.TEcalctype  = 'VW_ds'; % use the new TE estimator (Wibral, 2013)

% ACT estimation and constraints on allowed ACT(autocorelation time)
cfgTEP.trialselect = 'no';
cfgTEP.actthrvalue = 0;   % threshold for ACT
cfgTEP.minnrtrials = 5;   % minimum acceptable number of trials
cfgTEP.maxlag      = 1000;

% kernel-based TE estimation
cfgTEP.flagNei = 'Mass' ;           % neigbour analyse type
cfgTEP.sizeNei = 4;                 % neigbours to analyse

data_ = TEprepare(cfgTEP, data)

function [effsnrsite] =  create_semi_simulated_ECoG(fname,options)
% function to add a fake signal window to a SPM MEEG file.
%
% It takes as input the filename (full file name with path) of an epoched
% file in SPM MEEG format and options for the simulations
% options:
% -.isgood: whether to add only on 'good' channels (1) or all (0, default)
% -.times: vector of times in ms to define the onset and offset of
%           simulated signal (e.g. [50 200], default: all window)
% -.condA: conditions to add signal to, in cell array of strings 
%           (e.g. {'A'}, default: all)
% -.condB: condition to consider as 'noise', cell array of strings 
%           (e.g. {'B'}, default: baseline signal in [-minimum 0] of all conditions)
% -.newname: string with appended name for simulated files (default: 'Simulated')
% -.permute_chans: flag to permute channels randomly before adding signal
%           (default: 0)
% -.SNR: range of imposed SNR to simulate, each generating a separate file
% -.sparsity: range in % of sparsity, representing the % of channels to add
%           signal to (default: 5:5:100)
% -.window: substructure to specify the 'type' ('rect' or 'ramp'), the
%           parameters ('args', e.g. 3 for slope of ramp) and the smoothing
%           ('smooth.smooth': 1 or 0, 'smooth.win': value in ms) of the
%           window to add to the signal. Default is rectangular window, no
%           smoothing ('.window.type='rect'','.window.smooth.smooth=0')
%--------------------------------------------------------------------------
% Written by J. Schrouff, 26/04/2018, for the Laboratory of Behavioral and 
% Cognitive Neuroscience (LBCN, STanford University) and the Centre for
% Medical Image Computing (CMIC, University College London).
%
% Code released with no warranty or support. 
% This code was used to generate semi-simulated HFB power in ECoG
% recordings, in the 2 publications:
% (1)  Interpreting weight maps in terms of cognitive or clinical 
%      neuroscience: nonsense? J. Schrouff & J. Mourao-Miranda (submitted)
% (2) Face Information is Anatomically Localized but Distributed in Time 
%     Across the Human Temporal Cortex. J. Schrouff, S. Salehi, O. Raccah, 
%     V. Rangarajan, S. Baek, J. Mourão-Miranda, Z. Helili, A.L. Daitch &
%     J. Parvizi (submitted)
% Please cite (1) when using.


% Initialization: getting inputs
% -------------------------------------------------------------------------
% Select file and gather file name and path
if nargin<1 || isempty(fname)
    rawname = spm_select(1,'mat','Select file with epoched rest data');
else
    rawname = fname;
end
[gpath,name,ext] = spm_fileparts(rawname);

try
    D = spm_eeg_load(rawname);
catch
    error('create_semi_simulated_ECoG:nofile','No SPM MEEG file selected')
end

if nargin<2
    options = struct('isgood',0,...
        'times',[],...
        'condA',[],...
        'condB',[],...
        'newname','Simulated',...
        'permute_chans',0,...
        'SNR',2:0.5:10,...
        'sparsity',5:5:100);
end

% Load channel labels, either only good channels (isgood=1) or all (=0, default)
if ~isfield(options,'isgood') || isempty(options.isgood) || ~options.isgood
    label = chanlabels(D);
elseif options.isgood
    allchans = 1: nchannels(D);
    label = chanlabels(D,setdiff(allchans,D.badchannels));
end

% Choose time points to build simulated signal on, in ms (default: all)
if ~isfield(options,'times') || isempty(options.times)
    t1 = indsample(D,min(D.time));
    t2 = indsample(D,max(D.time));
else
    t1 = indsample(D,options.times(1)/1000);
    t2 = indsample(D,options.times(2)/1000);
end

% Choose 'simulated' conditions (condA)
if ~isfield(options,'condA') || isempty(options.condA)
    disp('Adding on all defined conditions')
    ia = indtrial(D,D.condlist,'good');
else
    ia = indtrial(D,options.condA,'good');
end

% Choose 'noise' conditions (condB): if empty, baseline in all conditions
% are used
if ~isfield(options,'condB') || isempty(options.condB)
    disp('Using baseline to define SNR')
    ib = indtrial(D,D.condlist,'good');
    tb1 = indsample(D,min(D.time));
    tb2 = indsample(D,0);
    if tb2<tb1 || isempty(tb2)
        error('Basline chosen but no time point before 0ms')
    end
else
    ib = indtrial(D,options.condB,'good');
    tb1 = t1;
    tb2 = t2;
end

% Name of simulated file
if isfield(options,'newname') && ~isempty(options.newname)
    newname = [options.newname,name,ext];
else
    newname = ['Simulated',name,ext];
end

% Permute channel labels randomly if specified
ig = indchannel(D,label);
if isfield(options,'permute_chans') && options.permute_chans
    sitesord = ig(randperm(length(ig)));
else
    sitesord = ig;
end

% Parameters to vary: SNR
if isfield(options,'SNR') && ~isempty(options.SNR)
    ratio = options.SNR;
else
    ratio = 2:0.5:10;
end

% Parameters to vary: Sparsity
if isfield(options,'sparsity') && ~isempty(options.sparsity)
    spars = round(options.sparsity*(length(ig)/100));
else
    spars = length(ig); % use all channels - i.e. 100% sparsity
end

% Parameters to vary: window
if ~isfield(options,'window') || isempty(options.window)
    options.window.type = 'rect';
    options.window.smooth.smooth = 0;
elseif isfield(options.window,'type') && strcmpi(options.window.type,'ramp') % Specify default slope for ramp signals
    if ~isfield(options.window,'args') || isempty(options.window.args)
        options.window.args = 3;
    end
end
if (isfield(options.window,'smooth') && isfield(options.window.smooth,'smooth')) && ...
        options.window.smooth.smooth %Smoothing required        
    if ~isfield(options.window.smooth,'win')
        options.window.smooth.smooth_win = 50;
    end    
else
    options.window.smooth.smooth = 0;
    options.window.smooth.win = 0;
end
smooth = options.window.smooth;
   

% Initialization: Computating std in 'noise' signal for SNR
% -------------------------------------------------------------------------
mest = zeros(length(ig),1);
for i = 1:length(ig)
    ci = mean(D(ig(i),tb1:tb2,ib),3);
    mest(i) = std(ci);
end


% Initialization: outputs
% -------------------------------------------------------------------------

% Input neural signal properties
insnrsite = zeros(length(ig),numel(spars),numel(ratio));
insparse = zeros(numel(spars),numel(ratio));
minSignal = zeros(length(ig),numel(spars),numel(ratio));

% Effective neural signal properties
effsnrsite = zeros(length(ig),numel(spars),numel(ratio));
meffSignal = zeros(length(ig),numel(spars),numel(ratio));


% Build simulated data sets
% -------------------------------------------------------------------------
for r=1:numel(ratio)
    % Make path for SNR value
    npath = [gpath,filesep,'Sparse_SNR_',num2str(ratio(r))];
    mkdir(gpath,['Sparse_SNR_',num2str(ratio(r))]);
    
    % For each sparsity ratio
    for rat = 1:length(spars)

        % Make path for sparsity combination
        npath2 = [npath,filesep,'Sparse_',num2str(round((spars(rat)/length(ig))*100))]; %Get folder name in effective % of sparsity
        mkdir(npath,['Sparse_',num2str(round((spars(rat)/length(ig))*100))]);
        cd(npath2)
        
        % Get raw data and clone
        D = spm_eeg_load(rawname);
        fnamedat =[npath2,filesep,newname];
        newD = clone(D,fnamedat);
        newD(:,:,:) = D(:,:,:);
        
        % Establish window parameters
        t = D.time;
        t_onset = t(t1);
        t_offset = t(t2); 
        insign = zeros(length(ig),length(t));
        for i=1:spars(rat)
            % Identify channel
            ic = ismember(ig,sitesord(i));
            
            % Build window 
            amp = mest(ic)*ratio(r); %Compute amplitude from SNR imposed            
            if strcmpi(options.window.type,'rect')
                [simul] = rectangular_window(t,t_onset,t_offset,amp);
            elseif strcmpi(options.window.type,'ramp')
                [simul] = ramp_window(t,t_onset,t_offset,options.window.args,amp,D.fsample);
            end
            
            % Smooth window if requested
            if smooth.smooth
                win_length = (smooth.win/1000)*D.fsample;
                gusWin= gausswin(win_length)/sum(gausswin(win_length));
                gusWin = gusWin';
                simul = conv(simul,gusWin,'same');
            end
            
            % Apply window to data
            news = zeros(size(D,2),numel(ia));
            for it = 1:numel(ia)
                news(:,it) = D(sitesord(i),:,ia(it))+ simul;
            end
            snrs(1,:,:) = news;
            newD(sitesord(i),:,ia) = snrs;
            
            % Estimate input neural signal properties
            insnrsite(ic,rat,r) = ratio(r);
            insign(ic,:) = simul;
            % Compute effective SNR (affected by window smoothness and shape)
            ci = mean(newD(sitesord(i),tb1:tb2,ib),3);
            mestd = std(ci);
            ca = mean(newD(sitesord(i),t1:t2,ia),3);
            ma = mean(ca);
            effsnrsite(ic,rat,r) =  ma/mestd;
        end
%         D = newD;
        save(newD);
%         D = spm_eeg_load(newD.fname);

        
        % Input sparsity pattern
        insparse(rat,r) = spars(rat)/length(ig);
        % Discrimination map - channel level
        insign = insign(:,t1:t2);
        minSignal(:,rat,r) = mean(insign,2);
        % Discrimination map - feature wise
        tor = mean(D(ig,t1:t2,ia),3) - mean(D(ig,t1:t2,ib),3);
        meffSignal(:,rat,r) = mean(tor,2);
         
%         % Plot resulting signal
%         fname = fullfile(npath2,newD.fname);
%         LBCN_plot_averaged_signal_epochs(fname,label);

        %Save temporary global results
        cd(gpath)
        save('Results_SNR_Sparsity.mat',...
            'effsnrsite',...
            'insnrsite',...
            'insparse',...
            'minSignal',...
            'meffSignal',...
            'sitesord')
        cd(npath)
    end
end

cd(gpath)




%--------------------------------------------------------------------------
% Subfunctions
%--------------------------------------------------------------------------
%Rectangular window
function [rect] = rectangular_window(t,t_onset,t_offset,amp)


sig1_inds = find(t<t_onset); % signal before onset
sig2_inds = find(t>=t_onset & t <= t_offset); %plateau
sig3_inds =find(t>t_offset); % signal after offset
sig1 = zeros(1,length(sig1_inds));
sig2 = amp*ones(1,length(sig2_inds));
sig3 = zeros(1,length(sig3_inds));

rect = [sig1,sig2,sig3];



%Ramp window
function [rect] = ramp_window(t,t_onset,t_offset,slope,amp,sampling_rate)

t_ramp = amp/slope;
sig1_inds = find(t<t_onset); % signal before onset
sig2_inds = find(t>=t_onset & t < t_onset+t_ramp); %ramp
sig3_inds = find(t>=t_onset+t_ramp & t<=t_offset); %plateau
sig4_inds =find(t>t_offset & t<=t_offset+t_ramp);
sig5_inds =find(t>t_offset+t_ramp);
sig1=zeros(1,length(sig1_inds));
sig2 = ((1:length(sig2_inds))/sampling_rate)*slope;
sig3 = amp*ones(1,length(sig3_inds));
sig4 = ((1:length(sig4_inds))/sampling_rate)*-slope +amp;
sig5 = zeros(1,length(sig5_inds));

rect = [sig1,sig2,sig3,sig4,sig5];
% --- set paths 
compName=lower(getenv('computername'));
switch compName
  case {'hh-i7'}
    sumResultDataPath='d:/hh/projects/ctx_propoSevo/data_IPSC/';
    dataPath='d:/_data/otc_ctx/SevoPropo/IPSC/';
    plotPath='d:/hh/projects/ctx_propoSevo/rawFig/';
  case {'hh64'}
    sumResultDataPath='d:/hh/projects/ctx_propoSevo/data_IPSC/';
    dataPath='e:/_data/otc_ctx/SevoPropo/IPSC/';
    plotPath='d:/hh/projects/ctx_propoSevo/rawFig/';
  otherwise
    error('unknown machine')
end

% --- graphics export
% graphics format within which to save plot (the usual input arguments into 
% matlab print function; set to [] for no plot)
gr.printas='-djpeg95';
gr.printas='-dpsc2';
gr.printas=[];
% directory to save figures in
gr.fDir=plotPath;

% --- layout of figure and subplots
% scaling factor for post processing, see help of labelscale.m
gr.scaleFac=.75;
% layout of figure (usual matlab choices 'tall', 'landscape', etc.)
gr.ornt='portrait';
% number of subplot rows (columns will be taken care of automatically)
gr.numSubPlotRows=5;
% font size 
gr.fontSz=8;
% default line width
gr.lineW=.25;

% --- appearance of continuous excerpt plot
% length of excerpt to be plotted (s)
gr.excLen=20;
% number of 'sweeps' per excerpt 
gr.nSweepPp=3;
% - offset between traces in case of several sweeps (original units of the
% recordings)
gr.dy=150;
% downsampling factor
gr.sampFac=2;

% --- appearance of overlaid cutouts plot
% number of cutouts per plot
gr.numCutout=20;
% peri-cutout interval (ms)
gr.cutoutIntv=[-5 50];
% interval for base line computation
gr.baselineIntv=[-3 -1];
% parameters and the limits within which the PSCs to plot must reside (**
% make sure you refer to fitted parameters)
gr.selectPSCPar={...
  'qFit', [.5 1.0];
  'amp', [50 150];
  };
% PSCs whose decay phase as given by gr.cutoutIntv contains a PSC with a
% *relative* amplitude larger than gr.disturbPSCNormAmpThresh will not be
% plotted
gr.disturbPSCNormAmpThresh=0.2;
% colors of averaged PSC cutout traces in overlay plot (will be cycled
% through)
gr.pscCol=[0 0 0; 0 0 1; .7 .7 .7];

%% ------------------------------------------------------------------------
% experimental date  
expDate='2016_03_02';
% raw data directory
ds.dDir=[dataPath expDate '_SET1/'];
% PSC results data file name
ds.sumResultFn=[sumResultDataPath '/IPSC_prop.mat'];
% raw data file list; column order:
% - file names WITHOUT EXTENSION 
% - start of the interval to plot (s)
% - title
ds.fList={...
  [expDate '_0014'], 92, 'control';...
  [expDate '_0017'], 100, '0.5 µM propofol';...  
  [expDate '_0018'], 10, 'washout'
  };
% channel; element order:
% - name
% - type
% - offset to be added to data (e.g. LJP) or 'baseline'
% - hi and lopass filter frequencies (nan to skip)
% - freq of line pickup filter (nan to skip)
% - clipping range
ds.chList={'IN 1', 'vc', 'baseline', [nan 2500], 50, [-4000 0]};
% seed value for random number generator (ensures that same cutouts are
% used each time PSCplotta is called)
ds.randSeed=3;

% call PSCplotta
PSCplotta(ds,gr)

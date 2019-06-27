dataP='e:\_data\otc_ctx\ACh\MikBFCtx\';
plotP='d:';

% --- line plot appearance 
% width of line plot representing data trace
gr.lineW=.5;
% font size (scale bar, title)
gr.fontSz=8;
% scaling factor for post processing, see help of labelscale.m
gr.scaleFac=1;
% marker size
gr.markSz=3;

% --- layout of figure and plot
% layout of figure (usual matlab choices 'tall', 'landscape', etc.)
gr.ornt='landscape';

% --- post-plot action settings
% graphics format within which to save plot (the usual input arguments into 
% matlab print function; set to [] for no plot)
% gr.printas='-dpsc2';
gr.printas='-djpeg97';
gr.printas=[];
% directory to save figures in
gr.fDir=[plotP '\_data\otc_ctx\ACh\MikBFCtx\hypblip\'];

% §§§ the complete set of independent parameters 
ds.indepPar=[0 1 2 -1];

% --- file-specific parameters that are common to the present set
% offset that should be added to data
ds.offset=-15;
% downsampling factor
ds.sampFac=5;
% pre-burst interval (ms) in which to average
% computed parameters and put out in periBuRC2
ds.preBuAvIntv=[-6000 0];

% *** a global variable that will accumulate results:
% - experiments in rows
% - independent pars in columns
% - each cell holds a 2D array with 
% -- as many rows as sweeps selected
% -- column 1: resistance (MOhm)
% -- column 2: tau (ms)
% -- column 3: pre-burst time (ms)
% -- column 4: post-burst time (ms)
% -- column 5: pre burst len (ms)
clear periBuRC periBuRC2;
global periBuRC periBuRC2
periBuRC=cell([0 numel(ds.indepPar)]);
% R in slice 1, tau in slice 2
periBuRC2=repmat(nan,[0 numel(ds.indepPar) 2]);


%% -- specific experiment -------------------------------
expDate='2011_11_09'; 
ds.dDir=[dataP expDate '_SET3\'];
% file list; column order:
% - file names WITHOUT EXTENSION 
% - independent parameter (e.g. drug concentration)
% - sweep selection (list of integers or 'a' for all in file)
% - the interval in POINTS on which analsis (e.g. fit) shall work
% - current amplitude injected (pA) 
fitStartPt=4135;
fitSpanPt=2400;
ds.fList={...
  [expDate '_0014'],  0, 1:50, [fitStartPt fitStartPt+fitSpanPt], -100;...
  [expDate '_0024'],  2, 1:50, [fitStartPt fitStartPt+fitSpanPt], -100;...
  [expDate '_0026'],  -1, 1:50, [fitStartPt fitStartPt+fitSpanPt], -100;...  
  };
% channel; column order:
% - name
% - type of analysis to perform
% - unused
% - unused
% - y limits in plot
ds.chList={'IN 0','hypblip', nan, nan, [-72 -40]};

% y limits of parameter plots
ds.coeff(1).ylim=[20 150]; % resistance
ds.coeff(1).tick=0:50:150;
ds.coeff(2).ylim=[1 12]; % time constant
ds.coeff(2).tick=0:2:14;

% call periburc
periburc(ds,gr)



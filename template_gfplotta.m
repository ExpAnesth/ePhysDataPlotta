% raw plot template using gfplotta
compName=lower(getenv('computername'));
switch compName
  case {'hh-i7'}
    dataDir = 'd:\_data\otc_ctx\ACh\ChR2BFCtx\';
    plotDir = 'd:\hh\projects\ctx_ACh\rawFig\';
  case {'hh64','hh-i5'}
    dataDir = 'e:\_data\otc_ctx\ACh\ChR2BFCtx\';
    plotDir = 'd:\hh\projects\ctx_ACh\rawFig\';
  case {'eval-lmb'}
    dataDir = 'i:\_data\otc_ctx\ACh\ChR2BFCtx\';
    plotDir = 'd:\hh\projects\ctx_ACh\rawFig\';
  otherwise
    error(['no paths yet set on machine ' compName]);
end

% --- line plot appearance 
% width of line plot representing data trace
gr.lineW=.25;
% font size (scale bar, title)
gr.fontSz=8;
% scaling factor for post processing, see help of labelscale.m
gr.scaleFac=1;
% color (order) of traces (in case more than one channel shall be plotted):
% 'k': all black
% a colormap (n by 3 array)
gr.colorOrder='k';
gr.colorOrder=coma('blackblue','n',2);

% --- layout of figure and plot
% layout of figure (usual matlab choices 'tall', 'landscape', etc
gr.ornt='landscape';
% number of 'sweeps' per excerpt (a value>1 means there will be 'line
% breaks')
gr.nSweepPp=3;

% --- post-plot action settings
% graphics format within which to save plot (the usual input arguments into 
% matlab print function; set to [] for no plot)
% gr.printas='-dpsc2';
gr.printas='-djpeg90';
gr.printas='-dpsc2';
gr.printas=[];
% directory to save figures in
gr.fDir=plotDir;

% downsampling factor
ds.sampFac=2;

%% ------------------------------------------------------------------------
% experimental date
expDate='2015_03_06';
% data directory
ds.dDir=[dataDir expDate '_SET1\'];
% file list; column order:
% - file names WITHOUT EXTENSION 
% - start of the interval to plot (s)
% - title
ds.fList={...
  [expDate '_0045'], 0,'';...
  [expDate '_0046'], 0,'';...
  [expDate '_0047'], 0,'';...  
  };
% channels; column order:
% - name
% - type
% - offset to be added to data (e.g. LJP) 
% - hi and lopass filter frequencies (nan to skip)
% - freq of line pickup filter (nan to skip)
% - clipping range
% - y 'ticks' (dotted lines)
ds.chList={...
%   'IN 8 LED', 'LED', 0, [nan nan], nan, [-inf inf], nan;...
  'IN 0', 'extra', 0, [nan nan], 50, [-inf inf], nan;...
  'IN 1', 'extra', 0, [nan nan], 50, [-inf inf], nan};
% length of excerpt to be plotted (s, same for all subplots)
gr.excLen=89.9;

gfplotta(ds,gr)

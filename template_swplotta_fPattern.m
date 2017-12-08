% --- set paths 
compName=lower(getenv('computername'));
switch compName
  case {'hh-i7'}
    dataPath='d:\_data\otc_ctx\ACh\';
    plotPath='d:\_data\otc_ctx\ACh\';
  case {'hh64'}
    dataPath='e:\_data\otc_ctx\ACh\';
    plotPath='e:\_data\otc_ctx\ACh\';
  case {'eval-lmb'}
    dataPath='h:\_data\otc_ctx\ACh\';
    plotPath='h:\_data\otc_ctx\ACh\';
  otherwise
    error('unknown machine')
end

% --- line plot appearance 
% line width in plots of data traces
gr.lineW=.5;
% font size 
gr.fontSz=8;
% scaling factor for post processing, see help of labelscale.m
gr.scaleFac=1;
% marker size
gr.markSz=3;

% --- layout of figure and plot
% layout of figure (usual matlab choices 'tall', 'landscape', etc.)
gr.ornt='tall';

% --- post-plot action settings
% graphics format within which to save plot (the usual input arguments into 
% matlab print function; set to [] for no plot)
% gr.printas='-dpsc2';
gr.printas='-djpeg95';
% gr.printas='-dpsc2';
gr.printas=[];

% directory to save figures in
gr.fDir=[plotPath '\X98\acute\stepCur_plot\'];

% file-specific parameters that are common to the present set
% offset that should be added to data
ds.offset=-17;
% downsampling factor
ds.sampFac=1;
% the minimal time (ms) the signal SLOPE has to be above threshold to be
% considered a spike
ds.minSpxSlopeWid=.1;

%% -- specific experiment -------------------------------
expDate='2012_02_16'; 
ds.dDir=[dataPath '\X98\acute\' expDate '_SET3\'];
% unique identifier of cell
ds.cellID='2012_02_16_s1c1';
% path & name of file into which computed parameters will be written via
% cellMaster
ds.cellTableFn=[plotPath '\X98\acute\cellTable_X98acute'];
% file list; column order:
% - file names WITHOUT EXTENSION 
% - sweep selection (list of integers or 'a' for all in file)
% - the interval in POINTS on which analysis (e.g. fit) shall work
% - start and step values of current injection (job 'depCurInj')
ds.fList={...
  [expDate '_0003'], 'a', [nan nan], [15 25];...
  [expDate '_0004'], 'a', [nan nan], [15 25];...  
  [expDate '_0005'], 'a', [nan nan], [15 25];...    
  };
% channel; column order:
% - name
% - type of analysis to perform
% - threshold for spike detection based on slope, mV/ms
% - spx analysis interval in ms 
% - y limits in plot (irrelevant for job 'depCurInj')
ds.chList={...
  'IN 1','depCurInj', 20, [0 399.7]+115.8 , [nan nan]};

% name of image file(s) (up to two) depicting neuron, cell with empty
% string otherwise
ds.imFile={''}

% call swplotta
swplotta(ds,gr)

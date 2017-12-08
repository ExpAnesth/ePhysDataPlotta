function PSCplotta(ds_in,gr_in)
% ** function PSCplotta(ds_in,gr_in) 
% is a high-level plotting function for post-synaptic currents (PSCs) which
% is designed to create i) raw data plots and ii) overlay plots of cutouts
% thereof with the following properties
% - one or more files per data set
% - one recording channel
% - equal axes in all subplots of the same sort
% The function is derived from gfplotta. In contrast to gfplotta, which
% produces 'working-style' kind of plots, PSCplotta is meant to generate
% publication-quality graphics and therefore features more manual control
% of the settings.
% For information on the input arguments ds_in and gr_in please see
% template_PSCplotta.m

% V 1.0 May 2017 Harald Hentschke

% to do: 
% - allow other than fitted parameters for selection of PSCs; make sure
% that either fitted or general parameters were selected

% ------------ PART 1: check input, adjust parameters, preparations
[ds,gr]=PSCplotta_defaultParams;
ds=checkFields(ds,ds_in);
gr=checkFields(gr,gr_in);  

% number of requested channels (counting duplicates, too)
nRequestChan=size(ds.chList,1);
% set rng to defined state
rng('default');
rng(ds.randSeed);

% --- graphics
labelscale('fontSz',gr.fontSz,'scaleFac',gr.scaleFac,'lineW',gr.lineW,'markSz',8); 
numFile=size(ds.fList,1);
if gr.numSubPlotRows<numFile
  warning('adapting number of subplot rows')
  gr.numSubPlotRows=numFile;
end
mkfig(1,'b'); clf
orient(gr.ornt);
figName=[ds.fList{1} '_PSCPlot' ];

% ------------ PART 2: load & preprocess data
% load PSC results file
load(ds.sumResultFn,'r');
% handles to raw and cutout subplots
rawSph=gobjects(1,numFile);
cutSph=gobjects(1,numFile);
for g=1:numFile
  fn=ds.fList{g,1};
  % retrieve information on raw data file
  [~,~,fi]=abfload([ds.dDir '\' fn '.abf'],'info');
  % ** convention: fi.si=original si, si=si after downsampling
  if gr.sampFac>1
    si=fi.si*gr.sampFac;
  else
    si=fi.si;
  end
  % --- i. select subset of detected/analyzed PSCs, if existent
  % locate current file in PSC results variables
  [resRow,resCol]=find(strcmp(fn,r.listExpFile));
  if ~isempty(resRow)
    existPSCResult=true;
    % number of PSCs (based on first parameter)
    nPSC=numel(r.pscr{resRow,resCol,strcmp(r.depPar,gr.selectPSCPar{1,1})});
    % logical index to PSCs fulfilling inclusion criteria
    OKIx=true(nPSC,1);
    for k=1:size(gr.selectPSCPar,1)
      ix=strcmp(r.depPar,gr.selectPSCPar{k,1});
      OKIx=OKIx & ...
      r.pscr{resRow,resCol,ix}>=gr.selectPSCPar{k,2}(1) & ...
        r.pscr{resRow,resCol,ix}<=gr.selectPSCPar{k,2}(2);
    end
    if any(OKIx) 
      % omit those PSCs from list followed by a sizeable PSC in their decay
      % time:
      % - index to tsls of detected PSCs
      ix1=strcmp(r.depPar,'allTsl');
      % - index to tsls of detected and fitted PSCs
      ix2=strcmp(r.depPar,'tsl');
      %       % - §§ make sure tsl is column array...
      %       r.pscr{resRow,resCol,ix2}=r.pscr{resRow,resCol,ix2}(:);
      % inter-PSC-intervals for all combinations of candidate PSCs and all
      % detected PSCs (** note new automatic arithmetic expansion here and
      % further below)
      tmpMat=r.pscr{resRow,resCol,ix1}'-r.pscr{resRow,resCol,ix2}(OKIx);
      % convert to logicals (true if within [0 gr.cutoutIntv(2)])
      tmpMat=tmpMat>0 & tmpMat<gr.cutoutIntv(2);
      % compute offending PSCs' amplitude relative to the candidate PSC's
      % amplitude
      ix1=strcmp(r.depPar,'allAmp');
      ix2=strcmp(r.depPar,'amp');
      tmpMat=tmpMat.*r.pscr{resRow,resCol,ix1}'./r.pscr{resRow,resCol,ix2}(OKIx);
      % eliminate entries in which the summed relative amplitude is above a
      % threshold of gr.disturbPSCNormAmpThresh
      OKIx=find(OKIx);
      OKIx(sum(tmpMat,2,'omitnan')>=gr.disturbPSCNormAmpThresh)=[];
    else
      OKIx=find(OKIx);
    end
    nOK=numel(OKIx);
    if nOK<gr.numCutout
      warning(['inclusion criteria for PSC cutouts are too strict, resulting in only ' int2str(nOK) ' cutouts']);
      if nOK<1
        warning('inclusion criteria for PSC cutouts are beyond values covered by data, picking regularly spaced tenth of PSCs')
        OKIx=unique(1:floor(nPSC/10):nPSC);
      end
    else
      % random selection, sorted
      OKIx=sort(OKIx(randperm(nOK,gr.numCutout)));
    end
    % extract tsl and other info needed further below
    psc.tsl=r.pscr{resRow,resCol,strcmp(r.depPar,'tsl')}(OKIx);
    psc.amp=r.pscr{resRow,resCol,strcmp(r.depPar,'amp')}(OKIx);
    % compute indexes to part of cutouts to be used for computation of base line
    psc.baselineIx=cont2discrete(gr.baselineIntv-gr.cutoutIntv(1),si/1000,'intv',true);
  else
    existPSCResult=false;
    warning('no PSC analysis available for current file');
    psc=[];
  end
  
  % --- ii. read and treat raw data
  % raw data excerpt limits in s
  intv=[ds.fList{g,2} ds.fList{g,2}+gr.excLen];
  if existPSCResult
    % load whole file because PSCs will be cut from it throughout its range
    d=abfload([ds.dDir '\' fn '.abf'],'channels',ds.chList(:,1)');
  else
    d=abfload([ds.dDir '\' fn '.abf'],'start',intv(1),'stop',intv(2),'channels',ds.chList(:,1)');
  end
  % ** abfload neither respects order of channels as specified in
  % ds.chList nor loads channels in duplicate, so let's do sth about it
  % - the order of channels as put out by abfload in d (no duplicates):
  outChanList=intersect(fi.recChNames,ds.chList(:,1),'stable');
  tmpIx=nan(1,nRequestChan);
  for ci=1:nRequestChan
    tmpIx(ci)=find(strcmp(ds.chList(ci,1),outChanList));
  end
  d=d(:,tmpIx);
  % lopass filter
  for ci=1:nRequestChan
    fifi=isfinite(ds.chList{ci,4});
    if fifi(2)
      d(:,ci)=lofi(d(:,ci),fi.si,ds.chList{ci,4}(2));
    end
  end  
  % downsample by sampFac, if requested...
  if gr.sampFac>1
    d=d(1:gr.sampFac:end,:);
    % si has already been adapted above
  end
  % finally, hipass filter
  for ci=1:numel(ds.chList(:,1))  
    fifi=isfinite(ds.chList{ci,4});
    if fifi(1)
      d(:,ci)=hifi(d(:,ci),si,ds.chList{ci,4}(1));
    end
  end
  % determine size
  [n1,n2]=size(d);
  offset=zeros(1,n2);
  for ci=1:n2
    % compute offset
    if isnumeric(ds.chList{ci,3})
      offset(ci)=ds.chList{ci,3};
    elseif ischar(ds.chList{ci,3}) && strcmp(ds.chList{ci,3},'baseline')
      offset(ci)= -median(d(:,ci));
    else
      warning('illegal value for offset')
    end
    % notch filter
    if isfinite(ds.chList{ci,5})
      d(:,ci)=elim_hum(d(:,ci),si,ds.chList{ci,5});
    end
  end
  
  % --- iii. cut out PSCs, bring in shape and cut d down to size
  if existPSCResult
    % produce cutouts
    psc.d=tsl2exc(d,si,{psc.tsl(:)},'win',gr.cutoutIntv);
    % subtract baseline (* note new automatic array expansion)
    psc.d=psc.d-mean(psc.d(psc.baselineIx(1):psc.baselineIx(2),:));
    % average unscaled
    psc.avD=mean(psc.d,2);
    % average of amplitude-scaled PSCs
    psc.avDScaled=mean(psc.d./psc.amp(:)',2);
    % time in ms
    psc.t=(0:size(psc.d,1)-1)'*(si/1000);
    % cut d down to size
    intvPts=cont2discrete(intv,si/1e6,'intv',1);
    d=d(intvPts(1):intvPts(2),:);
    n1=diff(intvPts)+1;
  end
  
  % --- iv. perform remaining operations on d
  % clip d, then add offset
  for ci=1:n2
    if isfinite(ds.chList{ci,6}(1))
      d(d(:,ci)<ds.chList{ci,6}(1),ci)=ds.chList{ci,6}(1);
    end
    if isfinite(ds.chList{ci,6}(2))
      d(d(:,ci)>ds.chList{ci,6}(2),ci)=ds.chList{ci,6}(2);
    end
    d(:,ci)=d(:,ci)+offset(ci);
  end
  % reshape & drop points at end
  tmpN1=gr.nSweepPp*floor(n1/gr.nSweepPp);
  if n1-tmpN1>0
    disp([int2str(n1-tmpN1) ' points had to be dropped at end']);
    d=d(1:tmpN1,:);
  end
  n1=tmpN1;
  d=reshape(d,n1/gr.nSweepPp,gr.nSweepPp);
  
  % ------------ PART 3: plot data -----------------------------------------
  % i. raw data
  rawSph(g)=subplot(gr.numSubPlotRows, 2,(g-1)*2+1);
  %   % pump up & reposition
  %   rexy('ax',gca,'xfac',1.1,'yfac',1.1);
  switch g
    case numFile
      % last plot: scale bar
      [~,~,~,ph]=pllplot(d,'si',si,'spacing','fixed','dy',gr.dy);
    otherwise
      % other plots: no scale bar
      [~,~,~,ph]=pllplot(d,'si',si,'spacing','fixed','dy',gr.dy,'noscb',1);
  end
  title(ds.fList{g,3},'interpreter','none','fontweight','normal');

  if existPSCResult
    % ii. psc cutouts
    curPSCCol=gr.pscCol(rem(g-1,size(gr.pscCol,1))+1,:);
    % - overlay, unscaled
    cutSph(g)=subplot(gr.numSubPlotRows, 4,(g-1)*4+3);
    ph=plot(psc.t,psc.d,'color',[.7 .7 .7]);
    hold on
    ph=plot(psc.t,psc.avD,'color','k','linewidth',2);
    axis tight off
    % overlay of averages of unscaled PSCs
    subplot(gr.numSubPlotRows,4,4);
    hold on
    plot(psc.t,psc.avD,'linewidth',2,'color',curPSCCol);
    % overlay of averages of scaled PSCs
    subplot(gr.numSubPlotRows,4,8);
    hold on
    plot(psc.t,psc.avDScaled,'linewidth',2,'color',curPSCCol);
  end
end
subpax(gcf,'spInd',rawSph);
cutSph(~isgraphics(cutSph))=[];
if ~isempty(cutSph)
  subpax(gcf,'spInd',cutSph);
  subplot(cutSph(end))
  utscaleb4('ms','pA','pos','lr');
  subplot(gr.numSubPlotRows,4,4);
  axis tight off
  utscaleb4('ms','pA','pos','lr');
  subplot(gr.numSubPlotRows,4,8);
  axis off
  niceyax
  utscaleb4('ms','pA','pos','lr');
end

drawnow
if ~isempty(gr.printas)
  print(gr.printas,'-r400','-painters',[gr.fDir figName]);
end


% ------------ LOCAL FUNCTIONS ------------------------------------------
% ------------ LOCAL FUNCTIONS ------------------------------------------
function f=checkFields(f,f_in)
s=fieldnames(f);
s_in=fieldnames(f_in);
sect=setdiff(s_in,s);
if ~isempty(sect)
  errordlg({'nonmatching field names(s): '; char(sect)});
end
% assign input values of input struct to current one
for g=1:numel(s_in)
  f.(s_in{g})=f_in.(s_in{g});
end

function [ds,gr]=PSCplotta_defaultParams
% --- graphics export
% graphics format within which to save plot (the usual input arguments into 
% matlab print function; set to [] for no plot)
gr.printas=[];
% directory to save figures in
gr.fDir='';

% --- layout of figure and subplots
% scaling factor for post processing, see help of labelscale.m
gr.scaleFac=1;
% layout of figure (usual matlab choices 'tall', 'landscape', etc.)
gr.ornt='tall';
% number of subplot rows (columns will be taken care of automatically)
gr.numSubPlotRows=3;
% font size 
gr.fontSz=8;

% --- appearance of continuous excerpt plot
% length of excerpt to be plotted (s)
gr.excLen=10;
% number of 'sweeps' per excerpt 
gr.nSweepPp=3;
% - offset between traces in case of several sweeps (original units of the
% recordings)
gr.dy=150;
% downsampling factor
gr.sampFac=1;
% line width
gr.lineW=.25;

% --- appearance of overlaid cutouts plot
% number of cutouts per plot
gr.numCutout=20;
% peri-cutout interval (ms)
gr.cutoutIntv=[-5 50];
% interval for base line computation
gr.baselineIntv=[-3 -1];
% parameters and the limits within which the PSCs to plot must reside 
gr.selectPSCPar={...
  'allAmp', [nan nan];
  };
% PSCs whose decay phase as given by gr.cutoutIntv contains a PSC with a
% *relative* amplitude larger than gr.disturbPSCNormAmpThresh will not be
% plotted
gr.disturbPSCNormAmpThresh=0.2;
% color of averaged PSC cutout trace
gr.pscCol='k';

% --- file specifics
% data directory
ds.dDir='';
% PSC results file name
ds.sumResultFn='';
% file list; column order:
% - file names WITHOUT EXTENSION 
% - start of the interval to plot (s)
% - title
ds.fList={...
  'file name', 0, '';...
  };
% channel; element order:
% - name
% - type
% - offset to be added to data (e.g. LJP) or 'baseline'
% - hi and lopass filter frequencies (nan to skip)
% - freq of line pickup filter (nan to skip)
% - clipping range
ds.chList={'channel name', 'signal type', 0, [nan nan], nan, [-inf inf]};
% seed value for random number generator (ensures that same cutouts are
% used each time PSCplotta is called)
ds.randSeed=0;


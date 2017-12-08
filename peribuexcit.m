function resColOrder=peribuexcit(ds_in,gr_in)
% ** function resColOrder=peribuexcit(ds_in,gr_in)
% a high-level plotting & analysis function which computes peri-burst
% parameters of cell excitability based on repetitive injection of
% sinusoidal ramp current injection (=sweep-based data).
% Details:
% - one or more files
% - ONE channel only 
% - one job: analyze results of sine+ramp current injection
% - computed parameters:
% -- number of spx
% -- noise and mean value of Em preceding stimulus

% *************************************************************************
% This function is a member of a family of functions which perform limited,
% specialized tasks on whole data sets (experiments, each of them composed
% of several files). In contrast to e.g. tslbatch and perievdeal, which are
% called once to work on a whole list of experiments, the functions are
% called for each experiment and place the results (if any) in global
% variables for further use, with varying degrees of safeguards against
% potential problems (duplicate calls on same data, etc.). See the specific
% functions. For orientation, the rough order of evolution is
% [gfplotta,swplotta] -> periburc -> peribuexcit -> pscdeal -> stimdeal
% *************************************************************************

% to do

global periBuEXCIT periBuEXCIT2 LISTEXP

% ------------ PART 0: define defaults
% --- description of global results variables:
% - periBuEXCIT is a cell array containing sweepwise information; each row
% is one experiment, treatments are in columns; order of columns within
% each cell:
resColOrder={'spx_stim','spx_nonstim','t_pre','t_post','t_in','E_M','bgAct'};
nResCol=numel(resColOrder);
% periBuEXCIT2 is a 3D array containing scalar information averaged or
% aggregated over all sweeps; each row is one experiment, treatments are in
% columns; order of SLICES:
% 1. mean (pre-burst)
% 2. mean (all)

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
gr.ornt='tall';

% --- post-plot action settings
% graphics format in which to save plot (the usual input arguments into 
% matlab print function; set to [] for no plot)
% gr.printas='-dpsc2';
gr.printas=[];
% directory to save figures in
gr.fDir=[''];

% --- file specifics
% data directory
ds.dDir=[''];
% file list; column order:
% - file names WITHOUT EXTENSION 
% - independent parameter (e.g. drug concentration)
% - sweep selection (list of integers or 'a' for all in file)
ds.fList={...
  'file name', 0, 'a';...
  };
% channel; column order:
% - name
% - type of analysis to perform
% - spx analysis interval in ms 
% - Em noise and mean analysis interval in ms 
% - y axis limits spx count
% - y axis limits Em percentiles
% - y limits in raw plot
ds.chList={'channel name','analysis type',[nan nan],[nan nan],[nan nan],[nan nan],[nan nan]};
% downsampling factor
ds.sampFac=1;
% offset that should be added to data
ds.offset=0;
% pre-burst interval (ms) in which to average computed parameters and put
% out in periBuEXCIT
ds.preBuAvIntv=[];
% the complete set of independent parameters
ds.indepPar=nan;

% ------------ PART 1: check input, adjust parameters, preparations
ds=checkFields(ds,ds_in);
gr=checkFields(gr,gr_in);  

job=ds.chList{1,2};

% make sure that independent parameters match
if ~isempty(setdiff([ds.fList{:,2}],ds.indepPar))
  error('at least one of the independent parameters specified for the current set of recordings is not in the global set (ds.indepPar)');
end

% --- graphics
labelscale('fontSz',gr.fontSz,'scaleFac',gr.scaleFac,'lineW',gr.lineW,'markSz',gr.markSz); 
nFile=size(ds.fList,1);
% main figure 
fh1=figure(1); clf, orient(gr.ornt);
% secondary figure
fh2=figure(2); clf, orient(gr.ornt);

scs=get(0,'screensize');
marg=round(scs(4)/40);
switch gr.ornt
  case 'landscape'
    set(fh1,'position',[scs(1)+marg  floor(scs(4)/2)-marg  scs(3)-2*marg  floor(scs(4)/2)-2*marg]);
    % set(fh2,'position',[scs(1)+marg  scs(2)+marg           scs(3)-2*marg  floor(scs(4)/2)-2*marg]);
  case 'tall'
    set(fh1,'position',[scs(1)+marg           scs(2)+2*marg  floor(scs(3)/2)-marg  scs(4)-5*marg]);
end
% figure 2 will be tall no matter what
set(fh2,'position',[floor(scs(3)/2)+marg  scs(2)+2*marg  floor(scs(3)/2)-marg  scs(4)-5*marg]);

% figure name: first of file names in set & channel name(s)
tmp=[ds.chList{:,1}];
figName=[ds.fList{1,1} '_' tmp(~isspace(tmp))];
% patch color
patchCol=[.8 .8 .8];

etslconst;

% ------------ PART 2: load, process and plot data
for g=1:nFile
  % clear anything that may have remained from previous file
  clear global bu evt 
  % index into column of global results variables periBuEXCIT and
  % periBuEXCIT2
  colIx=find(ds.indepPar==ds.fList{g,2});
  % ***********************************************************************
  if size(ds.chList,1)>1
    error('code can currently deal only with one channel');
  end
  chIx=1;
  % ***********************************************************************
  % channel name
  if g==1
    deblChName=ds.chList{1,1};
    deblChName=deblChName(~isspace(deblChName));
    % name of current experiment: first file + channel name
    curExpName=[ds.fList{1,1} '_' deblChName];
    % row index for current experiment into global results variable: if
    % experiment name exists already, point to corresponding row
    rowIx=find(strcmp(curExpName,LISTEXP));
    if isempty(rowIx)
      rowIx=size(periBuEXCIT,1)+1;
      LISTEXP{rowIx,1}=curExpName;
    elseif numel(rowIx)>1
      error('internal:duplicate experiments in LISTEXP');
    else
      disp(['overwriting values for ' curExpName]);
    end
    % ** this preallocates for new entries or wipes current ones
    periBuEXCIT2(rowIx,:,:)=nan;
    [periBuEXCIT(rowIx,:,:)]=deal({[nan]});
  end
  
  % determine number of sweeps and load as many as possible
  sweepList=ds.fList{g,3};
  [d,si,fi]=abfload([ds.dDir ds.fList{g,1} '.abf'],'info');
  disp(['**** number of available sweeps: ' int2str(numel(fi.sweepStartInPts))]);
  if isnumeric(sweepList)
    if ~isempty(setdiff(sweepList,1:numel(fi.sweepStartInPts)))
      warning('at least some of the requested sweeps do not exist - adjusting list')
      sweepList=intersect(sweepList,1:numel(fi.sweepStartInPts));
      % if no overlap between requested and existent sweeps sth is really
      % foul, so give up
      if isempty(sweepList)
        error('no overlap between requested and existing sweeps');
      end
    end
  else
    sweepList=1:numel(fi.sweepStartInPts);
  end
  % load raw data
  
%   % §§§ this is a hack: load all channels
%   [d,si]=abfload([ds.dDir ds.fList{g,1} '.abf'],'sweeps',sweepList);
  
  [d,si]=abfload([ds.dDir ds.fList{g,1} '.abf'],'channels',ds.chList(chIx,1),'sweeps',sweepList);
  
  
  % downsample, add offset
  d=d(1:ds.sampFac:end,:,:)+ds.offset;
  si=si*ds.sampFac;
  [n1 nCh nSw]=size(d);
  % *** the general assumption here is that there is no gap between
  % episodes, so check for this ***
  tmp=unique(diff(fi.sweepStartInPts));
  if ~(numel(tmp)==1 && floor(tmp/ds.sampFac)==n1)
    warning('nonzero gaps between episodes - some results/plots will be incorrect');
  end
  
  switch job
    case 'sineramp'
      % --------- preps 
      stimIntv=ds.chList{chIx,3};
      bgIntv=ds.chList{chIx,4};
      spxCountYLim=ds.chList{chIx,5};
      emPercYLim=ds.chList{chIx,6};
      excYLim=ds.chList{chIx,7};
      rawPlotSweepInd=unique(round(linspace(1,nSw,5)));
      % --------- load spx and bu data
      % load spx tsl
      load([ds.dDir ds.fList{g,1} '_' deblChName '_SPX_res.mat'],'evt');
      % rename
      tsl=evt.tsl{1};
      tslSweepIx=evt.sweepIx{1};
      % convert ts to points
      tslPts=cont2discrete(tsl,si/1000,'intv',0);
      % next, see whether bursts have been detected and saved for current file;
      % if so, load
      matFn=([ds.dDir ds.fList{g,1} '_' deblChName '_MP_res.mat']);
      isBuFileFound=exist(matFn,'file');
      if isBuFileFound
        load(matFn,'bu');
      end
      % --------- count spx
      % pick sweep-related spx, separate them into those during stimulus
      % and those not (could be before or after)
      spxCount=zeros(nSw,2);
      for h=1:nSw
        tmpTsl=tsl(tslSweepIx==sweepList(h));
        % inside
        spxCount(h,1)=numel(find(tmpTsl>stimIntv(1) & tmpTsl<=stimIntv(2)));
        % outside
        spxCount(h,2)=numel(tmpTsl)-spxCount(h,1);
      end
      % ---------------- alignment to bursts 
      % transfer spx results to cell and preallocate other columns with
      % nans
      periBuEXCIT{rowIx,colIx}=[spxCount nan([nSw nResCol-2])];
      if isBuFileFound 
        if ~isempty(bu.etsl) || ~isempty(bu.silentEtsl)
          tArr=discrete2cont(floor((repmat(fi.sweepStartInPts(sweepList)+ds.sampFac,1,2))/ds.sampFac),...
            si/1000)+repmat(stimIntv,nSw,1);
          % ** compute peri-burst subintervals
          tArr=etslsubdiv(bu.etsl,bu.silentEtsl,tArr);
          % implant time information in periBuEXCIT
          periBuEXCIT{rowIx,colIx}(1:nSw,3:5)=tArr(:,[1 2 3]);
        end
      end
      % in the rare cases of missing or empty etsl set pre-burst time to -1
      % so that these values make it into the averaging process
      if all(isnan(periBuEXCIT{rowIx,colIx}(:,3)))
        periBuEXCIT{rowIx,colIx}(:,3)=-1;
      end
      
      % --------- deal with raw data
      % reshape
      d=permute(d,[1 3 2]);
      bgIntvPts=cont2discrete(bgIntv,si/1000,'intv',1);
      % *** eliminate spx before determining percentiles of Em
      delIx=find(tslPts>=bgIntvPts(1) & tslPts<bgIntvPts(2));
      if ~isempty(delIx)
        % the indexes of affected sweeps
        tmpSweepIx=tslSweepIx(delIx);
        % unique values thereof
        uTmpSweepIx=unique(tmpSweepIx);
        % those currently loaded up
        uTmpSweepIx=intersect(uTmpSweepIx,sweepList);
        for sIx=1:numel(uTmpSweepIx)
          % index to tsl of current sweep
          tmpTslIx=delIx(tmpSweepIx==uTmpSweepIx(sIx));
          d(:,sweepList==uTmpSweepIx(sIx),1)=...
            etslexcsubst(d(:,sweepList==uTmpSweepIx(sIx),1),si,tsl(tmpTslIx),[-.5 4],[-1 4.5]);
        end
      end
      % percentiles of Em
      prctEm=prctile(d(bgIntvPts(1):bgIntvPts(end),:,1),[2.5 97.5 50]);
      % median Em
      periBuEXCIT{rowIx,colIx}(:,6)=prctEm(3,:)';

      if 1 
        % THE STANDARD
        % degree of ongoing activity (currently, percentile range)
        bgAct=diff(prctEm(1:2,:));

      elseif 1
        % hipass filter to extract fast signals
        tmpD=hifi(d(:,:,1),si,5);
        bgAct=diff(prctile(tmpD(bgIntvPts(1):bgIntvPts(end),:),[5 95]));
        
      elseif 0
        % §§§ hack #2: determine bgAct from properly filtered LFP data
        % (will only work if hack #1, loading all data, is in effect)
        d(:,:,2)=lofi(d(:,:,2),si,30);
        bgAct=diff(prctile(d(bgIntvPts(1):bgIntvPts(end),:,2),[5 95]));
      end
      
      periBuEXCIT{rowIx,colIx}(:,7)=bgAct';
      
      % --------------------- plot ----------------------------------------
      warning('off','MATLAB:Axes:NegativeDataInLogAxis');
      figure(fh1),
      % i. raw data (selected sweeps)
      rawSph(g)=subplot(3,nFile,g);
      dPlot=d(:,rawPlotSweepInd);
      dPlot(end:end+100,:)=nan;
      pllplot(dPlot(:),'si',si,'noscb',g<nFile);
      axis on
      grid on
      set(gca,'ylim',excYLim,'xtick',[]);
      title([ds.fList{g,1} ', ' ds.chList{chIx,1} ', ' num2str(ds.fList{g,2})],'interpreter','none');

      % ii. EM percentiles & spx counts
      subplot(3,nFile,nFile+g)
      % Em percentiles
      paH=patch([1:nSw nSw:-1:1]',[prctEm(1,:) fliplr(prctEm(2,:))]',patchCol);
      set(paH,'edgecolor',patchCol);
      axis([1 nSw emPercYLim]);
      ylabel('E_M');
      % ** second axis for spx counts
      ax2=axes('position',get(gca,'position'),'color','none',...
        'xaxisloc','bottom','yaxisloc','right','ygrid','on');
      ylabel('number of spx');
      hold on
      % in current version plot only stim-induced spx
      ph=plot(spxCount(:,1),'ko-');
      set(ph(1),'markerfacecolor','k');
      % mark sweeps to be displayed as raw data
      tmpIx=rawPlotSweepInd;
      ph=plot(tmpIx,spxCount(tmpIx,1),'mo');
      set(ph,'markersize',10);
      axis([1 nSw spxCountYLim])
      
      % iii. number of evoked spx depending on timing of stim to bursts
      ph=subplot(3,nFile,g+2*nFile); hold on
      % index to stimuli in desired pre-burst interval
      ix=periBuEXCIT{rowIx,colIx}(:,3)>ds.preBuAvIntv(1) & periBuEXCIT{rowIx,colIx}(:,3)<ds.preBuAvIntv(2);
      % add a little noise to the discrete spx counts to be able to see
      % different pre-burst episodes
      ns=rand(size(ix))*.1;
      plot(periBuEXCIT{rowIx,colIx}(:,3),periBuEXCIT{rowIx,colIx}(:,1)+ns,'ko-');
      plot(periBuEXCIT{rowIx,colIx}(:,4),periBuEXCIT{rowIx,colIx}(:,1)+ns,'ko-');
      % within-burst ones in red
      ph=plot(periBuEXCIT{rowIx,colIx}(:,5),periBuEXCIT{rowIx,colIx}(:,1)+ns,'ko-');
      set(ph,'color',[.8 .1 .1]);
      axis tight
      set(gca,'ylim',spxCountYLim,'ygrid','on');
      % line at t=0
      lh=line([0 0],spxCountYLim,'linestyle','--','color','r','linewidth',1.5);
      mn=nanmean(periBuEXCIT{rowIx,colIx}(ix,1));
      % plot this value as line
      lh=line(ds.preBuAvIntv,mn*[1 1],'color','g','linewidth',3);
      % ...and store in first slice
      periBuEXCIT2(rowIx,colIx,1)=mn;
      % second slice=spx averaged over all sweeps
      periBuEXCIT2(rowIx,colIx,2)=nanmean(periBuEXCIT{rowIx,colIx}(:,1));
      if g==1
        ylabel('number of spx');
        title([ds.fList{g,1} ', ' ds.chList{1,1} ', ' num2str(ds.fList{g,2})],'interpreter','none');
      end
      xlabel('peri-burst time (ms)');
      
      figure(fh2),
      
      % spx raster plot
      subplot(4,nFile,g);
      ph=plot(tslSweepIx,tsl,'k.');
      set(ph,'markersize',10);
      axis([sweepList(1) sweepList(end) ds.chList{chIx,3}])
      set(gca,'ygrid','on');
      
      % number of spx vs median Em
      corrSph1(g)=subplot(4,nFile,nFile+g);
      ph=plot(prctEm(3,:)',spxCount(:,1),'ko');
      set(ph,'markersize',3);
      nicexyax;
      xlabel('median E_M (mV)')
      ylabel('number of spx');
      % correlation coeff of this 
      [r,p,rlo,rup]=corrcoef(spxCount(:,1),prctEm(3,:)'); 
      smarttext(['r=' num2str(r(1,2),'%5.2f') '[' num2str(rlo(1,2),'%5.2f') ' ' num2str(rup(1,2),'%5.2f') ']' ]);
      
      % number of spx vs noise/bgAct
      corrSph2(g)=subplot(4,nFile,2*nFile+g);
      ph=plot(bgAct',spxCount(:,1),'ko');
      set(ph,'markersize',3);
      nicexyax;
      xlabel('noise in E_M (mV)')
      ylabel('number of spx');
      % correlation coeff of this 
      [r,p,rlo,rup]=corrcoef(spxCount(:,1),bgAct'); 
      smarttext(['r=' num2str(r(1,2),'%5.2f') '[' num2str(rlo(1,2),'%5.2f') ' ' num2str(rup(1,2),'%5.2f') ']' ]);
      
      % grand finale: number of spx, coded as the size of the symbol, in 2D
      % phase space defined by noise/bgAct of Em (abscissa) and Em (ordinate)
      corrSph3(g)=subplot(4,nFile,3*nFile+g);
      [ms,ar]=mapSpxCount2area(spxCount(:,1));
      col=repmat(.5,numel(ms),3);
      col(:,1)=min(ms/10,1);
      col(:,3)=1-min(ms/10,1);
      scH=scatter(bgAct',prctEm(3,:)',ar,col);
      nicexyax;
      grid on
      set(gca,'xscale','log');
      warning('off','MATLAB:Axes:NegativeDataInLogAxis');
      xlabel('noise in E_M (mV)')
      ylabel('E_M (mV)')
  end
end

warning('off','MATLAB:Axes:NegativeDataInLogAxis');
subpax(gcf,'spInd',corrSph1);
warning('off','MATLAB:Axes:NegativeDataInLogAxis');
subpax(gcf,'spInd',corrSph2);
warning('off','MATLAB:Axes:NegativeDataInLogAxis');
subpax(gcf,'spInd',corrSph3);

warning('off','MATLAB:Axes:NegativeDataInLogAxis');
if ~isempty(gr.printas),
  print(gr.printas,'-r400','-f1',[gr.fDir figName '_' job]);
  print(gr.printas,'-r400','-f2',[gr.fDir figName '_' job '2']);  
end
drawnow;
warning('on','MATLAB:Axes:NegativeDataInLogAxis');
      


% ------------ LOCAL FUNCTIONS ------------------------------------------
% ------------ LOCAL FUNCTIONS ------------------------------------------
function f=checkFields(f,f_in)
s=fieldnames(f);
s_in=fieldnames(f_in);
sect=setdiff(s_in,s);
if ~isempty(sect)
  errordlg({'nonmatching field names(s): '; strvcat(sect)});
end
% assign input values of input struct to current one
for g=1:numel(s_in)
  f.(s_in{g})=f_in.(s_in{g});
end

function [ms,ar]=mapSpxCount2area(spxCount)
% map spx count to SIZE of markers to be plotted in points
% - zero should be markersize 2
offs=2;
% - ten spx, which is about the max to be expected, shall be size 12
scale=[10 12];
% do it:
ms=spxCount*(scale(2)-offs)/scale(1)+offs;
% don't forget to convert from markersize to area, assuming a circle
ar=pi*(ms/2).^2;
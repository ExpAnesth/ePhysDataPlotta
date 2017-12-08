function periburc(ds_in,gr_in)
% ** function periburc(ds_in,gr_in) 
% a high-level plotting & analysis function which computes peri-burst 
% membrane parameters based on repetitive hyperpolarizing current injection
% (=sweep-based data).
% Details:
% - one or more files
% - ONE channel only 
% - two jobs; 'hypblip' computes values for each sweep individually and
% then averages (or combines) them; 'hypblip_avg' averages the sweeps
% before computing parameters
% 
% - computed parameters:
% -- membrane resistance and time constant (current clamp
%    measurements with hyperpolarizing current injection)
% 
% ** WARNING: the function blindly APPENDS data to global results variables
% periBuRC periBuRC2, meaning that there is no safeguard against duplicate
% entries whatsoever - see pscdeal.m for a solution
% ** NOTE: since July 2013, the code can save R and tau as computed during
% CONTROL conditions in a matfile using the \beastybites\cellmaster
% functions. The idea behind cellmaster is to have a global list of neurons
% recorded under control conditions and their properties so that we can
% characterize them based on many electrophysiological parameters. 

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
% - *** implement detrending sweeps! We need a new interval for that
% - filemaster should write values from hypblip_avg
% - user choice of type of fit
% - linear slopes won't do the trick with many data in which Ih is strong -
% find alternative to model this current

global periBuRC periBuRC2

% ------------------------------------------------------------------------- 
%          PART 1a: define gr and ds defaults
% -------------------------------------------------------------------------
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
% unique identifier of cell (last four letters: s=slice, c=cell)
ds.cellID='yyyy_mm_dd_sxcx';
% path & name of file into which computed parameters will be written via
% cellMaster
ds.cellTableFn='';
% file list; column order:
% - file names WITHOUT EXTENSION 
% - sweep selection (list of integers or 'a' for all in file)
% - the interval in POINTS on which analsis (e.g. fit) shall work
ds.fList={...
  'file name', 'a', [nan nan];...
  };
% channel; column order:
% - name
% - type of analysis to perform
% - unused
% - unused
% - y limits in plot
ds.chList={'channel name','analysis type', nan, nan, [nan nan]};
% downsampling factor
ds.sampFac=1;
% offset that should be added to data
ds.offset=0;
ds.preBuAvIntv=[];
% y limits of parameter plots
ds.coeff(1).ylim=[nan nan]; % amplitude
ds.coeff(1).tick=0:5:50;
ds.coeff(2).ylim=[nan nan];
ds.coeff(2).tick=0:5:40;

% the complete set of independent parameters
ds.indepPar=nan;
% row index for current experiment into global results variable
rowIx=size(periBuRC,1)+1;

% ------------------------------------------------------------------------- 
% PART 1b: check input, adjust parameters, preparations
% -------------------------------------------------------------------------
ds=checkFields(ds,ds_in);
gr=checkFields(gr,gr_in);  
% **
job=ds.chList{1,2};

% make sure that independent parameters match
if ~isempty(setdiff([ds.fList{:,2}],ds.indepPar))
  error('at least one of the independent parameters specified for the current set of recordings is not in the global set (ds.indepPar)');
end

% set to true if fit shall include a linear slope
doFitWithSlope=false;
% --- fitting parameters (hyperpolarizing current injection)
if doFitWithSlope
  fType = fittype(...
    'a*exp(-t/tau) + b*t + c',...
    'independent',{'t'},...
    'coefficients',{'a','tau','b','c'});
  fOpt = fitoptions(...
    'method','NonlinearLeastSquares',...
    'MaxIter',100,...
    'upper',[70 50 10 0],...
    'lower',[0 1 0 -100],... % !!do not allow negative slope!!
    'TolX',.01,...
    'display', 'off'...
    );
else
  fType = fittype(...
    'a*exp(-t/tau) + c',...
    'independent',{'t'},...
    'coefficients',{'a','tau','c'});
  fOpt = fitoptions(...
    'method','NonlinearLeastSquares',...
    'MaxIter',100,...
    'upper',[70 80 0],...
    'lower',[0 1 -140],...
    'TolX',.01,...
    'display', 'off'...
    );
end

etslconst;

% ------------------------------------------------------------------------- 
%                              PART 1c: graphics
% -------------------------------------------------------------------------
labelscale('fontSz',gr.fontSz,'scaleFac',gr.scaleFac,'lineW',gr.lineW,'markSz',gr.markSz); 
nFile=size(ds.fList,1);
% main figure
fh1=figure(1); clf, orient(gr.ornt);
scs=get(0,'screensize');
marg=round(scs(4)/40);
switch gr.ornt
  case 'landscape'
    set(fh1,'position',[scs(1)+marg floor(scs(4)/2)-marg scs(3)-2*marg floor(scs(4)/2)-2*marg]);
  case 'tall'
    set(fh1,'position',[scs(1)+marg scs(2)+2*marg floor(scs(3)/2)-marg scs(4)-5*marg]);
end
% main part of figure name: first of file names in set & channel name(s)
tmp=[ds.chList{:,1}];
figName=[ds.fList{1,1} '_' tmp(~isspace(tmp))];
% number of columns in main plot
nCol=1;
colOrd=get(gca,'colororder');

if strcmp('hypblip',ds.chList(:,2))
  % second figure (for overlaid sweeps and histograms)
  fh2=figure(2); clf, orient tall
  % figure 2 always tall
  set(fh2,'position',[scs(1)+marg  scs(2)+2*marg  floor(scs(3)/3)-marg  scs(4)-5*marg]);
  % third figure (plot of Rcell versus peri-burst time)
  fh3=figure(3); clf, orient tall
  set(fh3,'position',[floor(scs(3)/3)+marg  scs(2)+2*marg  floor(scs(3)/2)-marg  scs(4)-5*marg]);
  % preparations for histograms of Em amplitude and time constant
  if isnumeric(ds.fList{1,3}) && numel(ds.fList{1,3})<=30
    binW=.5;
  else
    binW=.2;
  end
  binA=0:binW:ds.coeff(1).ylim(2);
  binTau=0:binW:ds.coeff(2).ylim(2);
  
  % containers for axis x limits and corresponding axis handles of figure 1
  % so that after everything is plotted x limits can be homogenized
  timeArray=[];
  axArray=[];
  % same for figure 2
  axArray2=[];
end

% ------------------------------------------------------------------------- 
%              PART 2: load, process and plot data
% -------------------------------------------------------------------------
% preallocate periBuRC2 with Nans
periBuRC2(rowIx,:,:)=nan;
% *** loop over files
for g=1:nFile
  % index into column of global results variable periBuRC
  colIx=find(ds.indepPar==ds.fList{g,2});
  % analysis interval within sweeps in pts
  anIntv_pts=ds.fList{g,4};
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
  [d,si]=abfload([ds.dDir ds.fList{g,1} '.abf'],'channels',ds.chList(:,1),'sweeps',sweepList);
  % permute
  d=permute(d,[1 3 2]);
  d=d(1:ds.sampFac:end,:)+ds.offset;
  si=si*ds.sampFac;
  [n1 n2]=size(d);
  % don't forget to adjust analysis interval to downsampling factor
  anIntv_pts=ceil(anIntv_pts/ds.sampFac);
  % x values for fits
  fitX=(0:diff(anIntv_pts))'*si/1000;
  % ** detrend excerpts
  % define two intervals from which to compute linear trend: i) 95 % of
  % points up to start of analysis interval, and ii) an equivalent portion
  % of points at the end of each trace. From each of these, the median
  % will be computed, and the line connecting these will be subtracted from
  % the cutouts
  tmp=round(.9*anIntv_pts(1));
  x1=1:tmp;
  y1=median(d(x1,:));
  % mean of abscissa values
  x1=mean(x1);
  % same for second interval at end of trace
  x2=n1-tmp+1:n1;
  y2=median(d(x2,:)); 
  x2=mean(x2);
  slope=(y2-y1)/(x2-x1);
  
  
%   % §§ no good like this because all raw plots are off
%   
%   % **subtract slope and baseline as defined by median of first interval
%   d=d-(1:n1)'*slope-repmat(y1,[n1 1]);
  
  
  
  % see whether bursts have been detected and saved for current file;
  % if so, load
  deblChName=ds.chList{:,1};
  deblChName=deblChName(~isspace(deblChName));
  matFn=([ds.dDir ds.fList{g,1} '_' deblChName '_MP_res.mat']);
  isBuFileFound=exist(matFn,'file');
  if isBuFileFound
    load(matFn)
  end
  % define starting values
  if doFitWithSlope
    % fit starting pars ('coefficients',{'a','tau','b','c'});)
    fStart=[20  15  0 mean(d(anIntv_pts(2),:,:))];
  else
    % fit starting pars ('coefficients',{'a','tau','c'});)
    fStart=[20  15  mean(d(anIntv_pts(2),:,:))];
  end
  set(fOpt,'Startpoint',fStart);
  
  switch job
    case 'hypblip'
      % ------------------ PART 0: preparations --------------------
      % *** the general assumption here is that there is no gap between
      % episodes, so check for this (but let's not worry about gaps of less
      % than a sample point as such values most likely reflect rounding
      % errors anyways)
      tmp=unique(diff(fi.sweepStartInPts));
      % if ~(numel(tmp)==1 && floor(tmp/ds.sampFac)==n1)
      if ~(numel(tmp)==1 && isalmost(tmp/ds.sampFac,n1,.1))
        warndlg('nonzero gaps between episodes - some results/plots will be incorrect');
      end
      % extract portions of sweeps to be analyzed and make a 2D array
      fitD=d(anIntv_pts(1):anIntv_pts(2),:,:);
      fitX=(0:size(fitD,1)-1)'*si/1000;
      % from fitD determine percentiles for y limits if none were specified
      if isempty(ds.chList{5}) || any(~isfinite(ds.chList{5}))
        ylim=prctile(fitD(:),[.5 99.5]);
      else
        ylim=ds.chList{5};
      end
      % to make sure we're fitting the correct excerpt extract portions of
      % sweeps to be analyzed plus some pts at end and beginning, if
      % possible, plot, and overplot fitD
      deltaT=200;
      plotIntv=[max(anIntv_pts(1)-deltaT,1)  min(anIntv_pts(2)+deltaT,n1)];
      plotD=d(plotIntv(1):plotIntv(2),:,:);
      % show data
      figure(fh2),
      subplot(nFile,2,(g-1)*2+1), cla, hold on
      plot((plotIntv(1):plotIntv(2))*ds.sampFac,plotD,'k');
      plot((anIntv_pts(1):anIntv_pts(2))*ds.sampFac,fitD,'r');
      axis tight
      set(gca,'ylim',ylim);
      title([ds.fList{g,1} ', ' ds.chList{1,1} ', ' num2str(ds.fList{g,2})],'interpreter','none');
      if g==nFile
        xlabel('pts');
      end
      % ------------------ PART 1: fit & plot time course -----------------
      % plot concatenated data so we can marvel at the traces while
      % the time-consuming fitting business runs in background
      figure(fh1)
      sph=subplot(nFile*2,nCol,(g-1)*2+1);
      axArray=cat(1,axArray,sph);
      rexy('ax',gca,'xfac',1.15,'yfac',1.0);
      % time in s as x axis
      time=(0:n1*n2-1)*(si/1e6);
      timeArray=cat(1,timeArray,time([1 end])');
      plot(time,d(:),'k');
      axis tight
      set(gca,'ylim',ylim);
      title([ds.fList{g,1} ', ' ds.chList{1,1} ', ' num2str(ds.fList{g,2})],'interpreter','none');
      drawnow;

      trcFit=[];
      gof=[];
      progressbar;
      for h=n2:-1:1
        % transfer fit values into array
        [tmpf,tmpgof]=fit(fitX,fitD(:,h),fType,fOpt);
        trcFit(h,:)=coeffvalues(tmpf);
        progressbar((n2-h+1)/n2);
      end
      % compute Rcell (MOhm) from amplitude
      trcFit(:,1)=trcFit(:,1)/abs(ds.fList{g,5})*1000;
      % *** §§ kick out values beyond y limits
      kickIx=trcFit(:,1)<ds.coeff(1).ylim(1) | trcFit(:,1)>ds.coeff(1).ylim(2);
      trcFit(kickIx,1)=nan;
      kickIx=trcFit(:,2)<ds.coeff(2).ylim(1) | trcFit(:,2)>ds.coeff(2).ylim(2);
      trcFit(kickIx,2)=nan;
      % plot parameter values vs time (assumption: zero gap between
      % episodes)
      blipTime=(1:n2)'*n1*(si/1e6);
      timeArray=cat(1,timeArray,blipTime([1 end]));
      subplot(nFile*2,nCol,g*2);
      rexy('ax',gca,'xfac',1.15,'yfac',1.0);
      % plot amplitude (=resistance) and tau
      [ax,h1,h2]=plotyy(blipTime,trcFit(:,1),blipTime,trcFit(:,2));
      axArray=cat(1,axArray,ax');
      set(h1,'marker','o');
      set(h2,'marker','v');
      % grid on only for one axis; otherwise it'll be too crowded
      set(ax(1),'ylim',ds.coeff(1).ylim,'ytick',ds.coeff(1).tick,'xlim',[0 time(end)],...
        'gridlinestyle',':','ygrid','on');
      set(ax(2),'ylim',ds.coeff(2).ylim,'ytick',ds.coeff(2).tick,'xlim',[0 time(end)],...
        'gridlinestyle',':','ygrid','off');
      set(get(ax(1),'Ylabel'),'String','R');
      set(get(ax(2),'Ylabel'),'String','tau');

      % ------------------ PART 2: histograms of parameters ---------------
      nA=histc(trcFit(:,1),binA);
      nT=histc(trcFit(:,2),binTau);
      figure(fh2),
      subplot(nFile,2,g*2), cla, hold on
      [ax,h1,h2]=plotyy(binA,nA,binTau,-nT,@bar);
      set(h1,'facecolor',colOrd(1,:),'edgecolor',colOrd(1,:),'BarWidth',1.0);
      set(h2,'facecolor',colOrd(2,:),'edgecolor',colOrd(2,:),'BarWidth',1.0);
      set(ax(1),'xlim',binA([1 end]),'xgrid','on');
      set(ax(2),'xlim',binTau([1 end]),'xaxisloc','top');
      axArray2=cat(1,axArray2,ax);
      if g==nFile
        % labels for fig 2
        figure(fh2)
        set(get(ax(1),'xlabel'),'string','R (MOhm)');
        set(get(ax(2),'xlabel'),'string','tau (ms)');
        % homogenize x limits in fig 1
        set(axArray,'xlim',[min(timeArray) max(timeArray)]);
        % % homogenize y limits (fig 2, histograms, for each parameter separately)
        set(axArray2,'yTickMode','auto');
        subpax(fh2,'spInd',axArray2(:,1));
        subpax(fh2,'spInd',axArray2(:,2));
      end
      % ------------------ PART 3: alignment to bursts ------------------
      % transfer fit results to cell and preallocate other columns with
      % nans
      periBuRC{rowIx,colIx}=[trcFit(:,[1 2]) repmat(nan,[n2 3])];
      if isBuFileFound 
        if ~isempty(bu.etsl) || ~isempty(bu.silentEtsl)
          % §§ to do:
          % - fit, smooth, kick out outliers
          % - incorporate bu len (possibly 3D plot, possibly glm)

          % start points of analysis intervals (ms)
          aiStartTsl_pts=floor(fi.sweepStartInPts(sweepList)/ds.sampFac)+anIntv_pts(1);
          % same for ends of analysis intervals
          aiEndTsl_pts=aiStartTsl_pts+diff(anIntv_pts);
          tArr=discrete2cont([aiStartTsl_pts aiEndTsl_pts],si/1000,'intv',1);
          % ** compute peri-burst subintervals
          tArr=etslsubdiv(bu.etsl,bu.silentEtsl,tArr);
          % implant time information in periBuRC
          periBuRC{rowIx,colIx}(1:n2,3:5)=tArr(:,[1 2 4]);
        end
      end
      % in the rare cases of missing or empty etsl set pre-burst time to -1
      % so that these values make it into the averaging process
      if all(isnan(periBuRC{rowIx,colIx}(:,3)))
        periBuRC{rowIx,colIx}(:,3)=-1;
      end
      
      figure(fh3), 
      ix=periBuRC{rowIx,colIx}(:,3)>ds.preBuAvIntv(1) & periBuRC{rowIx,colIx}(:,3)<ds.preBuAvIntv(2);
      for pIx=1:2
        subplot(nFile,2,(g-1)*2+pIx),
        hold on
        plot(periBuRC{rowIx,colIx}(:,3),periBuRC{rowIx,colIx}(:,pIx),'ko-');
        plot(periBuRC{rowIx,colIx}(:,4),periBuRC{rowIx,colIx}(:,pIx),'ko-');
        axis tight
        set(gca,'ylim',ds.coeff(pIx).ylim,'ytick',ds.coeff(pIx).tick,'ygrid','on');
        % line at t=0
        lh=line([0 0],ds.coeff(pIx).ylim,'linestyle','--','color','r','linewidth',1.5);
        mn=nanmedian(periBuRC{rowIx,colIx}(ix,pIx));
        % plot this value as line
        lh=line(ds.preBuAvIntv,mn*[1 1],'color','g','linewidth',3);
        periBuRC2(rowIx,colIx,pIx)=mn;
        if pIx==1
          ylabel('R (MOhm)');
          title([ds.fList{g,1} ', ' ds.chList{1,1} ', ' num2str(ds.fList{g,2})],'interpreter','none');
          % ** write to cellTable if indep par==0 (control condition)
          if ds.fList{g,2}==0
            cellmaster_writepar(ds.cellTableFn,ds.cellID,{'R'},{mn});
          end
        else
          ylabel('tau (ms)');
          % ** write to cellTable if indep par==0 (control condition)
          if ds.fList{g,2}==0
            cellmaster_writepar(ds.cellTableFn,ds.cellID,{'tau'},{mn});
          end
        end
        if g==nFile
          xlabel('peri-burst time (ms)');
        end
      end
      
    case 'hypblip_avg'      

      % ------------------ PART 1: ----------------- 
      % 
      
      % separate sweeps into those preceding bursts (=during relative
      % quiescence) and those within bursts:
      % should there be no sweeps corresponding to bursts and silent
      % periods (either because ~isBuFileFound or because of choice of
      % parameters), by default define all sweeps as belonging to silent
      % periods
      pop(1).name='silent';
      pop(1).ix=1:n2;
      pop(2).name='burst';
      pop(2).ix=[];
      if isBuFileFound 
        if ~isempty(bu.etsl) || ~isempty(bu.silentEtsl)
          % start points of analysis intervals (ms)
          aiStartTsl_pts=floor(fi.sweepStartInPts(sweepList)/ds.sampFac)+anIntv_pts(1);
          % same for ends of analysis intervals
          aiEndTsl_pts=aiStartTsl_pts+diff(anIntv_pts);
          tArr=discrete2cont([aiStartTsl_pts aiEndTsl_pts],si/1000,'intv',1);
          % ** compute peri-burst subintervals
          tArr=etslsubdiv(bu.etsl,bu.silentEtsl,tArr);
          % index to sweeps in silent period in given interval before burst
          pop(1).ix=find(tArr(:,1)>ds.preBuAvIntv(1) & tArr(:,1)<ds.preBuAvIntv(2));
          % index to sweeps in bursts
          pop(2).ix=find(tArr(:,3)>0);
        end
      end
      % average the pops, fit, plot
      for j=1:numel(pop)
        if ~isempty(pop(j).ix)
          pop(j).mn=mean(d(:,pop(j).ix),2);
          pop(j).std=std(d(:,pop(j).ix),0,2);
          fitD=pop(j).mn(anIntv_pts(1):anIntv_pts(2));
          [tmpf,tmpgof]=fit(fitX,fitD,fType,fOpt);
          trcFit=coeffvalues(tmpf);
          % compute Rcell (MOhm) from amplitude
          trcFit(1)=trcFit(1)/abs(ds.fList{g,5})*1000;
          % plot
          deltaT=200;
          plotIntv=[max(anIntv_pts(1)-deltaT,1)  min(anIntv_pts(2)+deltaT,n1)];
          plotD=pop(j).mn(plotIntv(1):plotIntv(2));
          plotVar=pop(j).std(plotIntv(1):plotIntv(2));
          % plot
          figure(fh1),
          subplot(nFile,2,(g-1)*2+j), cla, hold on
          set(gca,'defaultlinelinewidth',2);
          % mean & std
          [h1,hp]=boundedline((plotIntv(1):plotIntv(2))*ds.sampFac,plotD,plotVar,'k-');
          % plot((plotIntv(1):plotIntv(2))*ds.sampFac,plotD,'k');
          % overplot with part used for fitting
          plot((anIntv_pts(1):anIntv_pts(2))*ds.sampFac,fitD,'r');
          % plot fit
          plot((anIntv_pts(1):anIntv_pts(2))*ds.sampFac,tmpf(fitX),'m');
          if j==1
            niceyax(2);
            ylim=get(gca,'ylim');
          else
            axis tight
            % use mean of fitD as center
            ylim=mean(fitD)+[-1 1]*diff(ylim)/2;
            set(gca,'ylim',ylim);
          end
          smarttext(['R=' num2str(trcFit(1)) ' M\Omega, \tau=' num2str(trcFit(2)) ' ms'],.95,.95);
          title([ds.fList{g,1} ', ' ds.chList{1,1} ', ' num2str(ds.fList{g,2}) ', ' pop(j).name],...
            'interpreter','none');
          if g==nFile
            xlabel('pts');
          end
        end
      end
      
  end
end

if strcmp('hypblip',ds.chList(:,2))
  chi=get(fh3,'chi');
  % in the new matlab version uimenus and uitoolbars are children, too
  chi=chi(strcmp(get(chi,'type'),'axes'));
  subpax(fh3,'spInd',chi(1:2:end))
  subpax(fh3,'spInd',chi(2:2:end))
end

if ~isempty(gr.printas),
  % print figure 1 in specified format and save as fig file
  print(gr.printas,'-r400','-f1',[gr.fDir figName '_' job]);
  %   saveas(fh1,[gr.fDir figName '_' job],'fig');
  if exist('fh2','var') && ishandle(fh2)
    print(gr.printas,'-f2',[gr.fDir figName '_' job '2']);
  end
  if exist('fh3','var') && ishandle(fh3)
    print(gr.printas,'-r300','-f3',[gr.fDir figName '_' job '3']);
  end
end
      


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


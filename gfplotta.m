function gfplotta(ds_in,gr_in)
% ** function gfplotta(ds_in,gr_in) 
% is a high-level function which creates plots of continuous ('gap-free')
% raw time series with the following properties
% - one or more files
% - one or more channels (optionally in different colors)
% - possibly 'line breaks'
% - equal axes in all subplots
% For information on the input arguments ds_in and gr_in please see
% template_gfplotta.m

% note: the 'signal type' characteristic plays only a role in the automatic
% scaling of signals (done via pllplot)

% ------------ PART 0: define defaults
% --- line plot appearance 
% width of line plot representing data trace
gr.lineW=.25;
% font size (title, scale bar, labels)
gr.fontSz=8;
% scaling factor for post processing, see help of labelscale.m
gr.scaleFac=1;
% color (order) of traces, two choices: must be either 'k' (all black) or a
% colormap (n by 3 array)
gr.colorOrder='k';

% --- layout of figure and plot
% layout of figure (usual matlab choices 'tall', 'landscape', etc.)
gr.ornt='landscape';
% length of excerpt to be plotted (s, same for all subplots)
gr.excLen=1;
% number of 'sweeps' per excerpt (a value>1 means there will be 'line
% breaks')
gr.nSweepPp=1;

% --- post-plot action settings
% graphics format in which to save plot (the usual input arguments into 
% matlab print function; set to [] for no plot)
% gr.printas='-dpsc2';
gr.printas=[];
% directory to save figures in
gr.fDir='';

% --- file specifics
% data directory
ds.dDir='';
% file list; column order:
% - file names WITHOUT EXTENSION 
% - start of the interval to plot (s)
% - title
ds.fList={...
  'file name', 0, ' ';...
  };
% channels; column order:
% - name
% - type
% - offset to be added to data (e.g. LJP) 
% - hi and lopass filter frequencies (nan to skip)
% - freq of line pickup filter (nan to skip)
% - clipping range
% - y 'ticks' (dotted lines)
ds.chList={'channel name', 'signal type', 0, [nan nan], nan, [-inf inf], nan;...
  };
% downsampling factor
ds.sampFac=1;

% ------------ PART 1: check input, adjust parameters, preparations
ds=checkFields(ds,ds_in);
gr=checkFields(gr,gr_in);  

% take into consideration earlier versions of gfplotta 
% - which did not require explicit title for each subplot
if size(ds.fList,2)>2
  isUserTitle=true;
else
  isUserTitle=false;
end
% - in which the offset was added to all channels
if size(ds.chList,2)<7
  errordlg('ds.chList has too few entries, possibly because you are using an outdated format - please see template_gfplotta.m for the current format, and particularly check whether there is an entry for ''offset to be added to data''');
end
% number of requested channels (counting duplicates, too)
nRequestChan=size(ds.chList,1);
% --- graphics
labelscale('fontSz',gr.fontSz,'scaleFac',gr.scaleFac,'lineW',gr.lineW,'markSz',8); 
if isnumeric(gr.colorOrder)&& size(gr.colorOrder,1)<size(ds.chList,1)
  warning('size of colormap not adequate - resetting to black traces');
  gr.colorOrder='k';
end
nPlot=size(ds.fList,1);
figure(1), clf, orient(gr.ornt);
scs=get(0,'screensize');
marg=round(scs(4)/40);
switch gr.ornt
  case 'landscape'
    set(gcf,'position',[scs(1)+marg floor(scs(4)/2)-marg scs(3)-2*marg floor(scs(4)/2)-2*marg]);
  case 'tall'
    set(gcf,'position',[scs(1)+marg scs(2)+2*marg floor(scs(3)/2)-marg scs(4)-5*marg]);    
end
figName=[ds.fList{1} '_gf' ];

% ------------ PART 2: load & preprocess data
% handle to subplots
sph=nan(1,nPlot);
for g=1:nPlot
  fn=ds.fList{g,1};
  intv=[ds.fList{g,2} ds.fList{g,2}+gr.excLen];
  
  % §§ maybe include option in ds to explicitly specify file type
  % §§ automatic expansion of data not yet implemented for formats other
  % than abf
  
  %   if exist([ds.dDir '\' fn '.mat'],'file')
  %     % ** watch out - a specific data format is expected - see abfmerge2mat
  %     [d,si]=matDload([ds.dDir '\' fn '.mat'],'start',intv(1),'stop',intv(2),'channels',ds.chList(:,1)');
  %   elseif exist([ds.dDir '\' fn '.raw'],'file')
  if exist([ds.dDir '\' fn '.raw'],'file')
    rawFInfo=rawinfo('filename',[ds.dDir '\' fn '.raw'],'print','no');
    fi=rawfi2genfi(rawFInfo,'channels',ds.chList{1,1});
    si=fi.si;
    d=rawload([ds.dDir '\' fn '.raw'],ds.chList(:,1)',[intv(1) intv(2)]*1000,rawFInfo);
    % put into array and convert to mV
    d=cat(2,d{:})/1000;
  else
    [~,~,fi]=abfload([ds.dDir '\' fn '.abf'],'info');
    % number of requested points per channel
    ppChan=diff(cont2discrete(intv*1000,fi.si/1000,'intv',0));
    % if requested interval is not in data, fill with median values of each
    % channel
    if ppChan>fi.dataPtsPerChan
      [d,si]=abfload([ds.dDir '\' fn '.abf'],'start',intv(1),'stop','e','channels',ds.chList(:,1)');
      % attach 
      d(fi.dataPtsPerChan+1:ppChan,:)=repmat(median(d),[ppChan-fi.dataPtsPerChan 1]);
    else
      [d,si]=abfload([ds.dDir '\' fn '.abf'],'start',intv(1),'stop',intv(2),'channels',ds.chList(:,1)');
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

    % old code below works only if no channel duplicates asked for
    %     tmpIx=nan(nRequestChan,1);
    %     for ci=1:nRequestChan
    %       tmpIx(ci)=find(strcmp(ds.chList(ci,1),fi.recChNames));
    %     end
    %     [nix,tmpIx]=sort(tmpIx);
    %     d=d(:,tmpIx);
    
  end
  
  % don't ask...
  % d(any(isnan(d),2),:)=[];
  
  % lopass filter
  for ci=1:nRequestChan
    fifi=isfinite(ds.chList{ci,4});
    if fifi(2)
      d(:,ci)=lofi(d(:,ci),si,ds.chList{ci,4}(2));
    end
  end  
  % downsample by sampFac, if requested...
  if ds.sampFac>1
    d=d(1:ds.sampFac:end,:);
    si=si*ds.sampFac;
  end
  % hipass filter
  for ci=1:numel(ds.chList(:,1))  
    fifi=isfinite(ds.chList{ci,4});
    if fifi(1)
      d(:,ci)=hifi(d(:,ci),si,ds.chList{ci,4}(1));
    end
  end
  % determine size
  [n1,n2]=size(d);
  
  for ci=1:n2
    % add offset
    d(:,ci)=d(:,ci)+ds.chList{ci,3};
    % notch filter
    if isfinite(ds.chList{ci,5})
      d(:,ci)=elim_hum(d(:,ci),si,ds.chList{ci,5});
    end
  end
  
  % clip data
  for ci=1:n2
    if isfinite(ds.chList{ci,6}(1))
      d(d(:,ci)<ds.chList{ci,6}(1),ci)=ds.chList{ci,6}(1);
    end
    if isfinite(ds.chList{ci,6}(2))
      d(d(:,ci)>ds.chList{ci,6}(2),ci)=ds.chList{ci,6}(2);
    end
  end

  % in case there are several channels, determine dy and y scale factor on
  % basis of whole sweep from first recording
  if g==1
    if n2>1
      % [nix,dy,ysf]=pllplot(d,'si',si,'noplot',1,'yscale',ds.chList(:,2),'spacing','percentile');
      [~,dy,ysf]=pllplot(d,'si',si,'noplot',1,'yscale',ds.chList(:,2));
      % pick max dy and scale up a bit so sweeps stand out
      dy=repmat([dy max(dy)*1.5],1,gr.nSweepPp);
      dy(end)=[];
    else
      dy='auto';
      ysf=ds.chList(:,2)';
    end
  end
  % reshape & drop points at end
  tmpN1=gr.nSweepPp*floor(n1/gr.nSweepPp);
  if n1-tmpN1>0
    disp([int2str(n1-tmpN1) ' points had to be dropped at end']);
    d=d(1:tmpN1,:);
  end
  n1=tmpN1;
  % 3D: sweeps in columns, channels along slices
  d=reshape(d(1:n1,:),[floor(n1/gr.nSweepPp),gr.nSweepPp,n2]);
  % now back to 2D
  d=permute(d,[1 3 2]);
  d=reshape(d,[floor(n1/gr.nSweepPp),gr.nSweepPp*n2]);

  % y scale factor
  if g==1
    ysf=repmat(ysf,1,gr.nSweepPp);
  end
  
  % ------------ PART 3: plot data -----------------------------------------
  sph(g)=subplot(nPlot,1,g);
  % pump up & reposition
  rexy('ax',gca,'xfac',1.1,'yfac',1.1);
  switch g
    case 1
      if ischar(dy) && strcmp(dy,'auto')
        [~,dy,ysf,ph]=pllplot(d,'si',si,'noscb',nPlot>1);
      else
        [~,dy,~,ph]=pllplot(d,'si',si,'spacing','fixed','dy',dy,'yscale',ysf,'noscb',nPlot>1);
      end
    case nPlot
      % last plot: scale bar
      [~,dy,~,ph]=pllplot(d,'si',si,'spacing','fixed','dy',dy,'yscale',ysf);
    otherwise
      % other plots: no scale bar
      [~,dy,~,ph]=pllplot(d,'si',si,'spacing','fixed','dy',dy,'yscale',ysf,'noscb',1);
  end
  if isnumeric(gr.colorOrder)
    for ti=1:numel(ph)
      set(ph(ti),'color',gr.colorOrder(mod(ti-1,n2)+1,:));
    end
  end
  
  % if 'ticks' are requested convert their values to actual y axis values,
  % factoring in y scale factor and offset of traces
  % - count number of legal tick entries
  tickI=zeros(n2,1);
  for ci=1:n2
    tickI(ci)=sum(isfinite(ds.chList{ci,7}));
  end
  if any(tickI)
    % shift subplot a trifle to left
    spp=get(gca,'pos');
    set(gca,'pos',[spp(1)*.75 spp(2:4)]);
    % collect tick values in nan-padded array
    ticks=nan(n2,max(tickI));
    for ci=1:n2
      ticks(ci,1:tickI(ci))=ds.chList{ci,7};
    end
    % build vector of trace offsets and replicate by max number of ticks
    % per channel
    tOffs=repmat(cumsum([0 dy])',[1 max(tickI)]);
    % replicate ticks by number of sweeps
    allTicks=repmat(ticks,[gr.nSweepPp 1]);
    % transform
    transTicks=allTicks./repmat(ysf',[1 max(tickI)])+tOffs;
    % linearize transformed and real ticks and determine index to finite
    % values
    finIx=isfinite(transTicks);
    allTicks=allTicks(finIx);
    transTicks=transTicks(finIx);
    % ensure that both allTicks and transTicks are column arrays (which is
    % not the case following lines above if there is only one channel, one
    % sweep but several ticks
    allTicks=allTicks(:);
    transTicks=transTicks(:);
    nTick=numel(allTicks);
    xCo=get(gca,'xlim')'*ones(1,nTick);
    line(xCo,repmat(transTicks,[1 2])','linestyle',':','color',[.2 .2 .2]);
    % label on right side
    text(xCo(2,:)' + diff(xCo(:,1))/100,transTicks,num2str(allTicks))
  end
  % title: file name & user-provided title (if existent)
  if isUserTitle
    % th=smarttext([fn ', ' ds.fList{g,3}],.01,.95);
    th=title([fn ', ' ds.fList{g,3}]);
  else
    % th=smarttext(fn,.01,.95);    
    th=title(fn);
  end
  set(th,'interpreter','none');
end
subpax(gcf,'spInd',sph);
drawnow
if ~isempty(gr.printas)
  print(gr.printas,'-r400',[gr.fDir figName]);
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


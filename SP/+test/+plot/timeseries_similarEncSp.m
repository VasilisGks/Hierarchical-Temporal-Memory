function plotTimeseries_similarEncSP(ts,encOnly,spOnly,encspOverlapHistory, tstep)
% Plots the timeseries and annotates it with the currently most-similar encodings and spatial pooler outputs

global ovlG;

if tstep==1
elseif tstep==2
  figure(1); clf;
  % Timeseries plot
  subplot(211);
  plot(ts); hold on;
  encOnlyG= scatter(encOnly, ts(encOnly), 40, 'MarkerEdgeColor','k','LineWidth',0.1,...
              'MarkerFaceColor','g');
  encOnlyG.XDataSource='encOnly'; encOnlyG.YDataSource='ts(encOnly)';
  spOnlyG= scatter(spOnly, ts(spOnly), 40, 'MarkerEdgeColor','k','LineWidth',0.1, ...
              'MarkerFaceColor','r');
  spOnlyG.XDataSource='spOnly'; spOnlyG.YDataSource='ts(spOnly)';
  overlapG= scatter(encspOverlapHistory{1}, ts(encspOverlapHistory{1}), 30,...
              'MarkerEdgeColor','k','LineWidth',0.1, ...
              'MarkerFaceColor',[0.85 0.85 0]);
  overlapG.XDataSource='encspOverlapHistory{1}'; overlapG.YDataSource='ts(encspOverlapHistory{1})';
  lineG= plot([tstep,tstep], [min(ts),max(ts)],'color','r');
  lineG.XDataSource='[tstep,tstep]';
  hold off;
  title('Spatial Pooler mapping property evaluation');
  legend([encOnlyG,spOnlyG,overlapG], 'top encodings', 'top SP', 'overlap');
  xlabel('timestep');
  grid minor;
  
  % Overlap plot
  subplot(212);
  %ovlG= stem(overlapLengths./maxOverlap*100,'-o','LineWidth',0.1,'Color','b',...
  %           'MarkerSize',4,'MarkerEdgeColor','g','MarkerFaceColor','r');
  ovlG= animatedline('LineStyle','-.','Marker','o','MarkerSize',4,'MarkerEdgeColor','g','MarkerFaceColor','r');
  axis([0,length(ts),0,100]);
  grid minor; xlabel('timestep'); ylabel('overlap');
  title('Percentage of overlapping encoder and SP SDRs');
else
  overlapFraction= length(encspOverlapHistory{1}) ./ encspOverlapHistory{2};
  set(groot,'CurrentFigure',1);
  test.plot.local_refreshdata(gcf,'caller');
  addpoints(ovlG,tstep,overlapFraction*100);
  drawnow;
end


function [ spHistory ] = testSPfun()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Test spatial pooler
%clear; close all;
 
ncolumns= 2048; ninput= [16*25,5*25,3*25];%256;
encParams= struct(...
  'ninput',ninput, ...
  'aspectRatio_2dtopo',25./24, ...
  'w',[21,27,27] ...%'buckets',220 ...             % 1-bits w=n-buckets+1
  );
spParams= struct(...
  'ncolumns',ncolumns,...
  'aspectRatio_2dtopo',2,...
  'initInhibitionRadius',3,...
  'globalInhibition',true,...
  'potentialDensity',0.85,...
  'activityDensity',0.02,...
  'stimulusThreshold',1,...
  'synPermActiveInc',0.05,...
  'synPermInactiveDec',0.008,...
  'connectPermThreshold',0.1,...
  'permMax',255,...
  'maxBoost',1,...
  'minOverlapDutyCycle',0.001,...
  'dutyCyclePeriod',1000 ...
  );

[ts,N]= test.inputTimeseries();
encParams= encoder.initPowerDay(encParams,ts);
sp= algo.SpatialPooler(sum(ninput),ncolumns,encParams.aspectRatio_2dtopo,spParams);

encHistory= zeros(sum(ninput),N); spHistory= false(ncolumns,N);
encspOverlapHistory= cell(N,2);
ts_plot= cellfun(@(x) x{1},ts);

for timestep=1:N
  tic;
  % encode, spatially pool
  encOut= encoder.powerDay(ts(timestep), encParams); % The encoder has no memory
  spOut= sp.compute(encOut, false,timestep);
  encHistory(:,timestep)= encOut; spHistory(:,timestep)= spOut;   % keep history for display
  
  % get the 10 SDR most similar to the last from the entire history
  similarEnc= test.topSimilarSDR(encHistory(:,1:timestep),10);
  similarSp= test.topSimilarSDR(spHistory(:,1:timestep),10);
  % find their overlap and differences
  encOnly= setdiff(similarEnc{2},similarSp{2});
  spOnly= setdiff(similarSp{2},similarEnc{2});
  encspOverlapHistory(timestep,:)= {intersect(similarEnc{2},similarSp{2}), length(similarEnc{2})};
  
  % show timeseries annotated with most similar encodings and SP outputs
  test.plot.timeseries_similarEncSp(ts_plot,encOnly,spOnly,encspOverlapHistory(timestep,:),timestep);
  %test.plot.encSpOut(encHistory, spHistory,timestep);
  % control framerate
  timestepDuration= toc; pause(0.01-timestepDuration);
end

ovp= cellfun(@(x) length(x), encspOverlapHistory(:,1)) ./ cellfun(@(x) x, encspOverlapHistory(:,2));
fprintf('Mean SP performance: [%.2f,%.2f]\n',mean(ovp(100:200)),mean(ovp(201:300)));
if exist('ovphist','var'), ovphist= [ovphist;[mean(ovp(100:200)),mean(ovp(201:300))]];
else, ovphist= [mean(ovp(30:130)),mean(ovp(131:230))]; end


end


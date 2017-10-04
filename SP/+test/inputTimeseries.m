function [ts,N]= inputTimeseries()

%N=300; ts= sin(6.283/6*(0:N-1)') + 0.1*randn(N,1);
%
timespan= 121:552; % similar to Matt's example
%timespan= 121:121+3*24;
N= length(timespan);
load('prediction-hourly.mat');
ts= cell(N,1);
for i= timespan
  ts{i-timespan(1)+1}= { ...
    reccenterhourly{i,2}, ...
    hours(timeofday(reccenterhourly{i,1})), ...
    isweekend(reccenterhourly{i,1}) ...
  };
end
%}

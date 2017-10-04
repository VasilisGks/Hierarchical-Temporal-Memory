function [sdr,b]= simpleArithm(x, params)
% Simple arithmetic encoding
% - x: point in timeseries
% - params {struct}:
%   - .ninput
%   - .inMax
%   - .inMin
%   - .buckets

b= floor((x-params.inMin) .* params.buckets./params.range) +1;
sdr= zeros(params.ninput,1);
sdr(b:b+params.w-1)= 1;
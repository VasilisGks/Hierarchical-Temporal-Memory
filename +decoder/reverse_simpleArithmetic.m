function mostLikely= reverse_simpleArithmetic(bucketProbDist, algorithm,params)
% Reverses the simple arithmetic encoder, given the same parameters. Inputs a probability
% distribution across all buckets and collapses it to a single arithmetic value, representing
% the most likely estimation in the timeseries' domain (like a defuzzifier)
% Arguments:
% - bucketProbDist [nbucket]: discrete probability of each bucket [0-1]
% - algorithm: {'mean','mode','highmean'} 'highmean' is the mean of the highest-estimated values
% - params {struct}:
%   - .ninput
%   - .inMax
%   - .inMin
%   - .buckets

assert(size(bucketProbDist,2)==1);          % probdist must be a column vector
bucketL= params.range/(2*params.buckets);
bucketC= params.inMin+bucketL+ (0:params.buckets-1)*2*bucketL;  % vector of bucket centers
bucketProbDist= bucketProbDist + eps(1/length(bucketC));        % normalize
bucketProbDist= bucketProbDist ./ sum(bucketProbDist);          % normalize

if strcmp(algorithm,'mode')           % Only the likeliest prediction, if it stands out
  if max(bucketProbDist)-min(bucketProbDist) > 0.2
    [~,idx]= max(bucketProbDist);
    mostLikely= bucketC(idx);
  else
    mostLikely= bucketC*bucketProbDist;
  end
elseif strcmp(algorithm,'mean')       % Mean of the entire distribution
  mostLikely= bucketC*bucketProbDist;
elseif strcmp(algorithm,'highmean')   % Mean of the most likely predictions (like "mean", but filtering)
  highprob= bucketProbDist >= quantile(bucketProbDist,0.85);
  bucketProbDist(highprob)= bucketProbDist(highprob) ./ sum(bucketProbDist(highprob));
  mostLikely= bucketC(highprob) * bucketProbDist(highprob);
else
  error(['Algorithm "',algorithm,'" doesn''t exist']);
end

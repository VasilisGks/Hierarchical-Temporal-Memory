function params= initSimpleArithm(params,ts)

params.inMax= max(ts); params.inMin= min(ts);
params.range= params.inMax - params.inMin +10*eps(params.inMax);
if isfield(params, 'buckets')
  params.w= params.ninput - params.buckets +1;
else
  params.buckets= params.ninput - params.w +1;
end

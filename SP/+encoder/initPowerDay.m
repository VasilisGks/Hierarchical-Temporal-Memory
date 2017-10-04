function params= initPowerDay(params,ts)

% Power encoder params
power= cellfun(@(x) x{1}, ts);
params.power.ninput= params.ninput(1);
params.power.w= params.w(1);
params.power.inMax= max(power); params.power.inMin= min(power);
params.power.range= params.power.inMax - params.power.inMin +12*eps(params.power.inMax);
params.power.buckets= params.power.ninput - params.power.w +1;

% Time encoder params
time= cellfun(@(x) x{2}, ts);
params.time.ninput= params.ninput(2);
params.time.w= params.w(2);
params.time.inMax= max(time); params.time.inMin= min(time);
params.time.range= params.time.inMax - params.time.inMin +12*eps(params.power.inMax);
params.time.buckets= params.time.ninput - params.time.w +1;

% Weekend encoder params
wkend= cellfun(@(x) x{3}, ts);
params.wkend.ninput= params.ninput(3);
params.wkend.w= params.w(3);
params.wkend.inMax= max(wkend); params.wkend.inMin= min(wkend);
params.wkend.range= params.wkend.inMax - params.wkend.inMin +12*eps(params.power.inMax);
params.wkend.buckets= params.wkend.ninput - params.wkend.w +1;

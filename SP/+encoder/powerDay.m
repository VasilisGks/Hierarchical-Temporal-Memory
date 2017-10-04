function [sdr,bucket]= powerDay(x,params)

x= x{1};
power= x{1};
[sdrPower,bucket]= encoder.simpleArithm(power,params.power);
time= x{2};
sdrTime= encoder.simpleArithm(time,params.time);
wkend= x{3};
sdrWkend= encoder.simpleArithm(wkend,params.wkend);

sdr= [sdrPower;sdrTime;sdrWkend];
function [ colsWithNoise ] = addNoise( cols,noiseVal,sparseCols)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%size1=colNum*cellNum
colsWithNoise=cols;
numOfCols=round(noiseVal*sparseCols);
 for i=1:numOfCols
   rng shuffle   %seed with time
   j=randi([1 length(cols)]); 
   colsWithNoise(j)=not(cols(j)); %change state of a col
 end

end


function [ activeSt , predictSt , learnSt ] = resetState( cellNum ,colNum ,learnSt )
size0=cellNum*colNum;

predictSt=zeros(1,size0);
%segmentSynapses=zeros(size0,size0); 
activeSt=zeros(1,size0);
learnSt=zeros(1,size0);
end


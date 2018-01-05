function [ segmentSynapses,synapsesWeight ] = growSynapses( segmentSynapses,synapsesWeight,learnStPrev,indx  )

maxNewSynapsesCount=32;
initialPermanence=0.21;
maxSynapsesPerSegment=32; %40 (OLD:val according to python impl.NEW=255/4=63.5~(63) = 63
connectedPermanence=0.5;

 if sum(segmentSynapses(indx,:))<maxSynapsesPerSegment
    newSynapses=abs(maxNewSynapsesCount-sum(synapsesWeight(indx,:)~=0));
    candidates=find(learnStPrev);
    numOfCands=min(length(candidates),newSynapses);
    for j=1:numOfCands %grow synpases to a subset of previous active columns
    pos=randi([1 numOfCands]); %random chose one candidate
    n=candidates(pos);
     if synapsesWeight(indx,n)==0 %new synapse grow 
       synapsesWeight(indx,n)=connectedPermanence; %initialize weight
     end
    end
  %   synapsesWeight=synapsesWeight+(synapsesWeight~=0)*initialPermanence;
end

end

function [backTrackingList, listPos,activeSegments,activeSt2,predictedSt,learnSt,synapsesWeight,segmentSynapses,numOfBurstingCols,numOfActiveCols] = tmpMemory2(tmParams,activeCols,predCellsPrev,learnStPrev,segmentSynapses,prevSegmentSynapses,synapsesWeight,prevActiveSt,learning,learnSt,backTrackingList,listPos)
import java.util.ArrayList
%input: segment - 2d array containing synapse connection (TOTAL_CELLS*TOTAL_CELLS{}
%input: segmentActive- 1d vector with state of each segment -
%Initial - segmentActive{}
%input: synapseWeight{}size=(col*cells,maxSynapsesPerSegment) - 2d array  containing weight for each synapse of
%synapseWeight : square array , holding value of each synapse
%TODO : More than one segments in each cell---
%(rows=numOfSegments*MaxSegmentsPerCell
%Parameter Initialization 
%TODO : TM in a Class'
%fprintf(tmParams.connectedPermanence)
size1=tmParams.colNum*tmParams.cellNum;
synapsesWeight = full(synapsesWeight);
%-----BackTracking Params---------
%-----BackTracking Params---------
%maxSegmentPerCell=255;
%maxSynapsesPerSegment=255;
%+ threshold , conn_perm,init_perm,inc_perm,dec_perm
tmpArray=(ones(1,tmParams.cellNum));      %multiply with column of predCellsPrevz
tmpArray2=ones(tmParams.cellNum*tmParams.colNum,1);
activeSt=uint8(zeros(tmParams.cellNum,tmParams.colNum));     %active state for cells
predictedSt=uint8(zeros(tmParams.cellNum*tmParams.colNum,1));  %flat array containing predictive states of all cells
%First Phase
activMat=tmpArray*predCellsPrev;  %Activate Cells with predictPreviousState=true
numOfBurstingCols=length(find((gather(activMat)==0).*activeCols'));
numOfActiveCols=length(find(activeCols));
act=find(activeCols); %indexes that map to active columns
%NEW ADD%%%%%%%%%%%%%%%
activeStPrev=reshape(prevActiveSt,tmParams.cellNum,tmParams.colNum); %TODW : !!Check an ontws swsto..
inactiveCells=((activeCols==0)'.*activeStPrev); %flat active state array(inactive cols) ?????
inactiveCells2=reshape(inactiveCells,1,size1);
%is a Row . Guarantee that active cols cells and synapses dont evolve in this operation
%Cells from inactive cols only 
inactArray=(repmat(inactiveCells2,size1,1)); %2d -- with size: (tmParams.colNum*tmParams.cellNum ^2) . 
inactSynapses=inactArray.*(synapsesWeight>=tmParams.connectedPermanence).*((sum(synapsesWeight'))>=(tmParams.activationThreshold)); %TODO :!!To tmParams.minThreshold kalyptei??/ 
synapsesWeight=synapsesWeight-(inactSynapses*tmParams.predictedDecrement);
%cells from inactive cols,with weights > 0,else 0s
%%%%END NEW%%%%%%
for i=1:length(act)   
  flag=0;
  index=act(i);  
  if(activMat(index)>=1) %at least 1 cell of col predicted
    activeSt(:,index) = uint8(predCellsPrev(:,index).*(tmpArray'));%.*actvSynaps'; %activate cells with true predictive state & active segment
    learnSt=reshape(learnSt,tmParams.cellNum,tmParams.colNum);
    learnSt(:,index)=activeSt(:,index); 
  else %activate all cells of col with mask array (ones)
    flag=1;
    activeSt(:,index) = tmpArray'; %activate all cells of current column  
  end
  activeSt2=reshape(activeSt,1,tmParams.colNum*tmParams.cellNum);
  %-----Backtracking Code-----------starts here-----------------
  learnStBT=reshape(learnSt,1,tmParams.cellNum*tmParams.colNum); %Flat learnSt vector
  if (backTrackingList.size)<tmParams.backTrackingWindow
          backTrackingList.add(sparse(double(learnStBT)));
  else
      backTrackingList.set(listPos,sparse(double(learnStBT)));
      if listPos<(tmParams.backTrackingWindow-1)
        listPos=listPos+1;   
      else
        listPos=0;
      end    
  end
  %-----Backtracking Code--------ends here----------------------------
if flag~=1   %IF current column isn't bursting
  for k=1:tmParams.cellNum  
     indx=((index-1)*tmParams.cellNum)+k;
     if activeSt(k,index)==1  %if active cell -grow new synapses
     [segmentSynapses,synapsesWeight]=growSynapses(tmParams,segmentSynapses,synapsesWeight,learnStPrev,indx,backTrackingList);
     end
     %if active segment - pos/neg reinforce active/inactive synapses 
     if (sum(and(segmentSynapses(indx,:),reshape(activeSt,size1,1)))>=tmParams.activationThreshold)
      synapsesWeight(indx,:)=synapsesWeight(indx,:)+(tmParams.incrWeight*(synapsesWeight(indx,:)>=tmParams.connectedPermanence))-(tmParams.decrWeight*(synapsesWeight(indx,:)<(tmParams.connectedPermanence)));
     %end
     %REMOVED
      %Punish matching segments 
     elseif(sum(and(segmentSynapses(indx,:),prevActiveSt))>=tmParams.minThreshold)
      synapsesWeight(indx,:)=synapsesWeight(indx,:)-(tmParams.decrWeight*(synapsesWeight(indx,:)>=tmParams.connectedPermanence));   
     end  
     %%REMOVED
  end     
else   %IF columns bursted
 if (flag==1)
 id=(index-1)*tmParams.cellNum+1; %Now finding cell with most active synapses in busted columns
 %Multiply -Weights for each cell with mask Array(ones(..)) with -all
 %Positives weights
 bestMatchCell=(synapsesWeight(id:(id+(tmParams.cellNum-1)),:)>=0)*tmpArray2;%Mult Arrs:(cell*(size1)*(size1))
 %Find Best Matching cell
 [val,pos]=max(bestMatchCell);  %Cell in burst col with highest weight sum
 id2=(id-1)+pos; %best matching cell index 
 learnSt=reshape(learnSt,1,tmParams.cellNum*tmParams.colNum);
 learnSt(id2)=1; %TODO : Keep or remove?
 [segmentSynapses,synapsesWeight]=growSynapses(tmParams,segmentSynapses,synapsesWeight,learnStPrev,id2,backTrackingList);
 synapsesWeight(id2,:)=synapsesWeight(id2,:)+(tmParams.incrWeight*(synapsesWeight(id2,:)>=tmParams.connectedPermanence))-(tmParams.decrWeight*(synapsesWeight(id2,:)<tmParams.connectedPermanence));
 end  
end
end 
activeSt2=reshape(activeSt,1,tmParams.colNum*tmParams.cellNum);
%%%%%%
segmentSynapses=synapsesWeight>=tmParams.connectedPermanence; %update weight locally
actvSynapses=and(segmentSynapses,repmat(activeSt2,size1,1));
activeSegments=sum(actvSynapses');
%sum up active cells of synapses
%%END NEW
predictedSt=(activeSegments>=tmParams.activationThreshold); %NEWWW.*(activeSt2==0); %Enter pred.state to cells that are not active in current t
cellsWithActvSegment=find(activeSegments);
activeSt2=reshape(activeSt,1,tmParams.colNum*tmParams.cellNum);
predCellsPrevFlat=reshape(predCellsPrev,tmParams.colNum*tmParams.cellNum,1);
xorArray=xor(predCellsPrevFlat,activeSt2');
punishArray=find(xorArray==1);
for i=1:length(punishArray)
    %negative reinforce synapses for cells wrong predicted
  j=punishArray(i);
  synapsesWeight(j,:)=synapsesWeight(j,:)-((tmParams.decrWeight*((synapsesWeight(j,:)>=tmParams.connectedPermanence)).*(synapsesWeight(j,:)~=0)));%-((tmParams.decrWeight*(activeSt2==0)).*(synapsesWeight(j,:)~=0));
end
reinfArray=find(xorArray==0);
for i=1:length(reinfArray)
  j=reinfArray(i);            
  synapsesWeight(j,:)=synapsesWeight(j,:)+((tmParams.incrWeight*(synapsesWeight(j,:)>=tmParams.connectedPermanence)))-((tmParams.decrWeight*(synapsesWeight(j,:)<tmParams.connectedPermanence)).*(synapsesWeight(j,:)~=0));
end
%NEW According to latest TM ,BAMI .LINE 49-54 .Punish inactive synapses
%last step . All existing synapses reinforced by initalPermanence value
segmentSynapses=synapsesWeight>=tmParams.connectedPermanence; %activate synapses above threshold
synapsesWeight=synapsesWeight+(synapsesWeight~=0)*tmParams.initialPermanence;
synapsesWeight=((synapsesWeight<=1).*synapsesWeight) + (synapsesWeight>1); %weights limit value  Weight_Borders1
%update weights locally
actvSynapses=and(segmentSynapses,repmat(activeSt2,size1,1));
activeSegments=sum(actvSynapses');
%sum up active cells of synapses
predictedSt=(activeSegments>=tmParams.activationThreshold); %NEWWW.*(activeSt2==0); %Enter pred.state to cells that are not active in current t

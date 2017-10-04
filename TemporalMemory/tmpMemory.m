function [activeSegments,activeSt2,predictedSt,learnSt,synapsesWeight,segmentSynapses,numOfBurstingCols,numOfActiveCols] = tmpMemory( activeCols,predCellsPrev,learnStPrev,segmentSynapses,prevSegmentSynapses,synapsesWeight,prevActiveSt,learning,learnSt)
%input: segment - 2d array containing synapse connection (TOTAL_CELLS*TOTAL_CELLS{}
%input: segmentActive- 1d vector with state of each segment -
%Initial - segmentActive{}
%input: synapseWeight{}size=(col*cells,maxSynapsesPerSegment) - 2d array  containing weight for each synapse of
%synapseWeight : square array , holding value of each synapse
%TODO : More than one segments in each cell---
%(rows=numOfSegments*MaxSegmentsPerCell
%Parameter Initialization 
%TODO : TM in a Class
colNum = 512;         %total number of columns in SP
cellNum=8;           %cells per column 
size1=colNum*cellNum;
activationThreshold=6;
initialPermanence=0.21;
connectedPermanence=0.65;
minThreshold=2;
maxNewSynapsesCount=10;
predictedDecrement=0.002;    
incrWeight=0.07;             %reinforcement of the currently active segment
decrWeight=0.07; 
%maxSegmentPerCell=255;
%maxSynapsesPerSegment=255;
%+ threshold , conn_perm,init_perm,inc_perm,dec_perm
tmpArray=ones(1,cellNum);      %multiply with column of predCellsPrevz
tmpArray2=ones(cellNum*colNum,1);
activeSt=zeros(cellNum,colNum);     %active state for cells
predictedSt=zeros(cellNum*colNum,1);  %flat array containing predictive states of all cells

%First Phase
activMat=tmpArray*predCellsPrev;  %Activate Cells with predictPreviousState=true
numOfBurstingCols=length(find((activMat==0).*activeCols'));
numOfActiveCols=length(find(activeCols));

act=find(activeCols); %indexes that map to active columns
for i=1:length(act)   
  flag=0;
  index=act(i);  
  if(activMat(index)>=1) %at least 1 cell of col predicted
    activeSt(:,index) = predCellsPrev(:,index).*(tmpArray');%.*actvSynaps'; %activate cells with true predictive state & active segment
    learnSt=reshape(learnSt,cellNum,colNum);
    learnSt(:,index)=activeSt(:,index); 
  else %activate all cells of col with mask array (ones)
    flag=1;
    activeSt(:,index) = tmpArray'; %activate all cells of current column  
  end
  activeSt2=reshape(activeSt,1,colNum*cellNum);
if flag~=1   %IF current column isn't bursted
  for k=1:cellNum  
     indx=((index-1)*cellNum)+k;
     if activeSt(k,index)==1  %if active cell -grow new synapses
     [segmentSynapses,synapsesWeight]=growSynapses( segmentSynapses,synapsesWeight,learnStPrev,indx);
     end
     %if active segment - pos/neg reinforce active/inactive synapses 
     if (sum(and(segmentSynapses(indx,:),activeSt2))>=activationThreshold)
      synapsesWeight(indx,:)=synapsesWeight(indx,:)+(incrWeight*(segmentSynapses(indx,:)==1))-(decrWeight*(segmentSynapses(indx,:)<connectedPermanence));
     %Punish matching segments 
     elseif(sum(and(segmentSynapses(indx,:),activeSt2))>=minThreshold)
      synapsesWeight(indx,:)=synapsesWeight(indx,:)-(decrWeight*(synapsesWeight(indx,:)>0));   
     end  
  end     
else   %IF columns bursted
 if (flag==1)
 id=(index-1)*cellNum+1; %Now finding cell with most active synapses in busted columns
 %Multiply -Weights for each cell with mask Array(ones(..)) with -all
 %Positives weights
 bestMatchCell=(synapsesWeight(id:(id+(cellNum-1)),:)>0)*tmpArray2;%Mult Arrs:(cell*(size1)*(size1))
 %Find Best Matching cell
 [val,pos]=max(bestMatchCell);  %Cell in burst col with highest weight sum
 id2=(id-1)+pos; %best matching cell index 
 learnSt=reshape(learnSt,1,cellNum*colNum);
 learnSt(id2)=1; %TODO : Keep or remove?
 [segmentSynapses,synapsesWeight]=growSynapses( segmentSynapses,synapsesWeight,learnStPrev,id2);
 synapsesWeight(id2,:)=synapsesWeight(id2,:)+(incrWeight*(segmentSynapses(id2,:)==1))-(decrWeight*(segmentSynapses(id2,:)<connectedPermanence));
 end  
end
end 
activeSt2=reshape(activeSt,1,colNum*cellNum);
%%%%%%
segmentSynapses=synapsesWeight>=connectedPermanence; %update weight locally
actvSynapses=and(segmentSynapses,repmat(activeSt2,size1,1));
activeSegments=sum(actvSynapses');
%sum up active cells of synapses
%%END NEW
predictedSt=(activeSegments>=activationThreshold); %NEWWW.*(activeSt2==0); %Enter pred.state to cells that are not active in current t
cellsWithActvSegment=find(activeSegments);
activeSt2=reshape(activeSt,1,colNum*cellNum);
predCellsPrevFlat=reshape(predCellsPrev,colNum*cellNum,1);
xorArray=xor(predCellsPrevFlat,activeSt2');
punishArray=find(xorArray==1);
for i=1:length(punishArray)
    %    synapsesWeight(i,:)=synapsesWeight(i,:)+tmpArrayIncr;
    %negative reinforce synapses for cells wrong predicted
  j=punishArray(i);
  synapsesWeight(j,:)=synapsesWeight(j,:)-((decrWeight*((segmentSynapses(j,:)==1)).*(synapsesWeight(j,:)~=0)));%-((decrWeight*(activeSt2==0)).*(synapsesWeight(j,:)~=0));
end
reinfArray=find(xorArray==0);
for i=1:length(reinfArray)
  j=reinfArray(i);            
  synapsesWeight(j,:)=synapsesWeight(j,:)+((incrWeight*(segmentSynapses(j,:)==1)))-((decrWeight*(segmentSynapses(j,:)<connectedPermanence)).*(synapsesWeight(j,:)~=0));
end
%NEW According to latest TM ,BAMI .LINE 49-54 .Punish inactive synapses
%last step . All existing synapses reinforced by initalPermanence value
synapsesWeight=synapsesWeight+(synapsesWeight~=0)*initialPermanence;
synapsesWeight=((synapsesWeight<=1).*synapsesWeight) + (synapsesWeight>1); %weights limit value  Weight_Borders1
segmentSynapses=synapsesWeight>=connectedPermanence; %activate synapses above threshold
%update weights locally
actvSynapses=and(segmentSynapses,repmat(activeSt2,size1,1));
activeSegments=sum(actvSynapses');
%sum up active cells of synapses
predictedSt=(activeSegments>=activationThreshold); %NEWWW.*(activeSt2==0); %Enter pred.state to cells that are not active in current t
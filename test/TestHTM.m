import java.util.ArrayList
backTrackingList = ArrayList();
%Temporal Memory Parameters..
tmParams = struct(...
    'colNum',1200,...
    'cellNum',15,...
    'activationThreshold',12,...
    'minThreshold',11,...
    'connectedPermanence',0.5,...
    'initialPermanence',0.21,...
    'sparsity',0.02,...
    'maxNewSynapsesCount',28,...
    'maxSynapsesPerSegment',32,... 
    'incrWeight',0.11,...
    'decrWeight',0.11,...
    'predictedDecrement',0.02,...
    'backTrackingWindow',12 ...
    );
tic
[ spHistory] = testSPfun();
%with testSP call
%load('cols.mat');
[numofSPcols,N]=size(spHistory);
size1=tmParams.colNum*tmParams.cellNum;

totalNum=tmParams.cellNum*tmParams.colNum;
columns=zeros(tmParams.colNum,1);
sparsity=0.02;
sparseCols=round(tmParams.colNum*sparsity);
predCellsPrev=((zeros(tmParams.cellNum,tmParams.colNum)));
predCellsPrevFlat=int8((zeros(1,tmParams.cellNum*tmParams.colNum)));
learnSt=int8((zeros(1,tmParams.cellNum*tmParams.colNum)));
learnStPrev=((zeros(1,tmParams.cellNum*tmParams.colNum)));
prevActiveSt=((zeros(1,tmParams.cellNum*tmParams.colNum)));
segmentSynapses=gpuArray(uint8((zeros(totalNum,totalNum))));
synapsesWeight=single((zeros(totalNum,totalNum)));
prevSegmentSynapses={}; 
x=[];y=[];z=[];w=[];prec=[];rec=[];anomaly_score=[];MAPE=[];pred_vals=[];act_vals=[];accur=[];pred_err=[];MAPE=zeros(N,1);
ProbArray=[];MSError=[];
learning=1; %enable learning
activeSynapses=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
enableNoise=0;
%Decoder Paramaters
ninput= [120,6*25,3*25];%256;
encParams= struct(...
  'ninput',ninput, ...
  'aspectRatio_2dtopo',25./24, ...
  'w',[27,27,27] ...%'buckets',220 ...             % 1-bits w=n-buckets+1
  );
%Window for anomaly likelihood
w=500;w0=15;
mt=[];var=[];mtBar=[];L=[];  
%%
bucketNum=ninput(1)-encParams.w(1)+1; %number of buckets for power consumption
pred_array=zeros(bucketNum,N);
[ts,N]= test.inputTimeseries();
encParams= encoder.initPowerDay(encParams,ts);
decoderAlpha=0.06;prediction_timesteps=5;
dec= decoder.SDRClassifier(prediction_timesteps, decoderAlpha,tmParams.cellNum*tmParams.colNum, encParams.power.buckets);
%%%%%%%%%%%%%%%%
predError=0;
energyMat=load('prediction-hourly.mat');
for stepNum=1:N
    if stepNum==1
     listPos = 0;
    end
    columns=spHistory(1:tmParams.colNum,stepNum);
    %[activeSt2,predictedSt,learnSt ] = resetState( tmParams.cellNum,tmParams.colNum,learnSt );
    %Calling TM function
    [backTrackingList,listPos,activeSynapses,activeSt2,predictedSt,learnSt,synapsesWeight,segmentSynapses,numOfBurstingCols,numOfActiveCols ] = tmpMemory2(tmParams, columns,predCellsPrev,learnStPrev,segmentSynapses,prevSegmentSynapses,synapsesWeight,prevActiveSt,learning,learnSt,backTrackingList,listPos);
    learnStPrev=learnSt';
    [encOut,targetBucket]= encoder.powerDay(ts(stepNum), encParams);
    predictedSt=gather(predictedSt);
    pred= dec.predict(predictedSt,targetBucket, true);
    pred_array(:,stepNum)=pred;
    %[~,bucketIndx]=max(pred);
   % [~,targetBucket]=max(pred);
   x= decoder.reverse_simpleArithmetic(pred,'mode',encParams.power); 
   energyConsRange=max(energyMat.reccenterhourly.kw_energy_consumption)-min(energyMat.reccenterhourly.kw_energy_consumption);%+10*eps(max(energyMat.reccenterhourly.kw_energy_consumption));
   %Finding LMS for predictions   
   %384 = n - w + 1(ninput(1) - w(1) + 1)
  % x=targetBucket*energyConsRange/encParams.power.buckets+min(energyMat.reccenterhourly.kw_energy_consumption)-10*eps(max(energyMat.reccenterhourly.kw_energy_consumption));
   pred_vals(stepNum)=x;act_vals(stepNum)=energyMat.reccenterhourly.kw_energy_consumption(stepNum+125);
   MSError(stepNum)=abs((energyMat.reccenterhourly.kw_energy_consumption(stepNum+125)-x));
   prevActiveSt=double(activeSt2);
    if length(sum(predCellsPrev))>0
    [predError,accuracy,recall,precis]=findAccuracy(activeSt2,reshape(predCellsPrev,1,tmParams.cellNum*tmParams.colNum));
    end
    predCellsPrev=reshape(predictedSt,tmParams.cellNum,tmParams.colNum);
    %Store Prediction error
    %And Accuracy for predictions in TM
    x(stepNum)=stepNum;accur(stepNum)=accuracy;anomaly_score(stepNum)=(numOfBurstingCols/numOfActiveCols);
    z(stepNum)=stepNum;pred_err(stepNum)=predError;prec(stepNum)=precis;rec(stepNum)=recall;  
     a1=abs(act_vals-pred_vals); 
   MAPE(stepNum) = sum(a1)/sum(act_vals);

    %Measure of Anomaly Likelihood  
  if mod(stepNum,10)==0
	      fprintf(string('MAPE='))
        a1=abs(act_vals-pred_vals); sum(a1)/sum(act_vals)	
  end
  
  
   if stepNum>w-1
    t=stepNum-(0:(w-1));
    mt(stepNum)=sum(anomaly_score(t))/w;   %Mean
    var(stepNum)=sqrt(sum(power((anomaly_score(t)-mt(stepNum)),2))/(w-1) ); %Variance
    t0=stepNum-(0:(w0-1));
    mtBar(stepNum)=sum(anomaly_score(t0)/w0);
    L(stepNum)=1-qfunc((mtBar-mt) / var);
   end
   %%%%%%%%%%%
   

end

  save('/home/gkitsasv/hotgym/OLDSP/act_vals.mat','-v7.3');
  save('/home/gkitsasv/hotgym/OLDSP/pred_vals.mat','-v7.3');
  save('/home/gkitsasv/hotgym/OLDSP/MAPE.mat','-v7.3');
  save('/home/gkitsasv/hotgym/OLDSP/anomaly_score.mat','-v7.3');    

RMSerror=sqrt(sum(MSError.^2)/N) %Computing Root Mean Squared Error
  %End of Script
 time_elapsed = toc;
 fprintf('Time elapsed : %u' ,time_elapsed);
 MSE_tr=0.5*sum(MSError.^2)
a1=abs(act_vals-pred_vals); sum(a1)/sum(act_vals)
 

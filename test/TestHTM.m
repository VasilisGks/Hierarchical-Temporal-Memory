tic
[ spHistory] = testSPfun();
%with testSP call
%load('cols.mat');
[numofSPcols,N]=size(spHistory);
colNum = 950;%original total number of columns in SP        %total number of columns in SP
cellNum=11;            %cells per column 
size1=colNum*cellNum;
activationThreshold=5;
initialPermanence=0.21;
connectedPermanence=0.5;
minThreshold=3;
maxNewSynapsesCount=8; 
totalNum=cellNum*colNum;
columns=zeros(colNum,1);
sparsity=0.02;
sparseCols=round(colNum*sparsity);
predCellsPrev=((zeros(cellNum,colNum)));
predCellsPrevFlat=int8((zeros(1,cellNum*colNum)));
learnSt=int8((zeros(1,cellNum*colNum)));
learnStPrev=((zeros(1,cellNum*colNum)));
prevActiveSt=((zeros(1,cellNum*colNum)));
segmentSynapses=gpuArray(uint8((zeros(totalNum,totalNum))));
synapsesWeight=single((zeros(totalNum,totalNum)));
prevSegmentSynapses={}; 
x=[];y=[];z=[];w=[];prec=[];rec=[];anomaly_score=[];MAPE=[];pred_vals=[];act_vals=[];accur=[];pred_err=[];
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
  'w',[30,27,27] ...%'buckets',220 ...             % 1-bits w=n-buckets+1
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
dec= decoder.SDRClassifier(prediction_timesteps, decoderAlpha, cellNum*colNum, encParams.power.buckets);
%%%%%%%%%%%%%%%%
predError=0;
energyMat=load('prediction-hourly.mat');

%Backtracking Code goes here
for stepNum=1:N
    columns=spHistory(1:colNum,stepNum);
    [ activeSt2,predictedSt,learnSt ] = resetState( cellNum,colNum,learnSt );
    %Calling TM function
    [activeSynapses,activeSt2,predictedSt,learnSt,synapsesWeight,segmentSynapses,numOfBurstingCols,numOfActiveCols ] = tmpMemory( columns,predCellsPrev,learnStPrev,segmentSynapses,prevSegmentSynapses,synapsesWeight,prevActiveSt,learning,learnSt);
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
    [predError,accuracy,recall,precis]=findAccuracy(activeSt2,reshape(predCellsPrev,1,cellNum*colNum));
    end
    predCellsPrev=reshape(predictedSt,cellNum,colNum);
    %Store Prediction error
    %And Accuracy for predictions in TM
    x(stepNum)=stepNum;accur(stepNum)=accuracy;anomaly_score(stepNum)=(numOfBurstingCols/numOfActiveCols);
    z(stepNum)=stepNum;pred_err(stepNum)=predError;prec(stepNum)=precis;rec(stepNum)=recall;  
    MAPE=(act_vals(stepNum)-pred_vals(stepNum))/act_vals(stepNum);
    %Measure of Anomaly Likelihood
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
RMSerror=sqrt(sum(MSError.^2)/N) %Computing Root Mean Squared Error
  %End of Script
 time_elapsed = toc;
 fprintf('Time elapsed : %u' ,time_elapsed);
 MSE_tr=0.5*sum(MSError.^2)
a1=abs(act_vals-pred_vals); sum(a1)/sum(act_vals)
 

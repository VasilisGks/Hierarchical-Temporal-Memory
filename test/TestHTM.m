[ spHistory] = testSPfun();
%with testSP call
%load('cols.mat');
[numofSPcols,N]=size(spHistory);
colNum = 512;%original total number of columns in SP        %total number of columns in SP
cellNum=8;            %cells per column 
size1=colNum*cellNum;
activationThreshold=6;
initialPermanence=0.21;
connectedPermanence=0.65;
minThreshold=2;
maxNewSynapsesCount=7; 
totalNum=cellNum*colNum;
columns=zeros(colNum,1);
sparsity=0.02;
sparseCols=round(colNum*sparsity);
predCellsPrev=zeros(cellNum,colNum);
predCellsPrevFlat=zeros(1,cellNum*colNum);
learnSt=zeros(1,cellNum*colNum);
learnStPrev=zeros(1,cellNum*colNum);
prevActiveSt=zeros(1,cellNum*colNum);
segmentSynapses=zeros(totalNum,totalNum);
synapsesWeight=zeros(totalNum,totalNum);
prevSegmentSynapses={}; 
x=[];y=[];z=[];w=[];prec=[];rec=[];anomaly_score=[];pred_vals=[];act_vals=[];
ProbArray=[];pred_array=zeros(266,N);
learning=1; %enable learning
activeSynapses=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
enableNoise=0;
%Decoder Paramaters
ninput= [12*25,3*25,3*25];%256;
encParams= struct(...
  'ninput',ninput, ...
  'aspectRatio_2dtopo',25./18, ...
  'w',[35,40,35] ...%'buckets',220 ...             % 1-bits w=n-buckets+1
  );
[ts,N]= test.inputTimeseries();
encParams= encoder.initPowerDay(encParams,ts);
decoderAlpha=0.01;prediction_timesteps=2;
dec= decoder.SDRClassifier(prediction_timesteps, decoderAlpha, cellNum*colNum, encParams.power.buckets);
%%%%%%%%%%%%%%%%
predError=0;
energyMat=load('prediction-hourly.mat');
for stepNum=1:N
    columns=spHistory(1:colNum,stepNum);
    [ activeSt2,predictedSt,learnSt ] = resetState( cellNum,colNum,learnSt );
    %Calling TM function
    [activeSynapses,activeSt2,predictedSt,learnSt,synapsesWeight,segmentSynapses,numOfBurstingCols,numOfActiveCols ] = tmpMemory( columns,predCellsPrev,learnStPrev,segmentSynapses,prevSegmentSynapses,synapsesWeight,prevActiveSt,learning,learnSt);
    learnStPrev=learnSt';
    [encOut,targetBucket]= encoder.powerDay(ts(stepNum), encParams);
    pred= dec.predict(predictedSt,targetBucket, true);
    pred_array(:,stepNum)=pred;
    %[~,bucketIndx]=max(pred);
    [~,targetBucket]=max(pred);
    energyConsRange=max(energyMat.reccenterhourly.kw_energy_consumption)-min(energyMat.reccenterhourly.kw_energy_consumption)+10*eps(max(energyMat.reccenterhourly.kw_energy_consumption));
   %Finding LMS for predictions   
   x=targetBucket*energyConsRange/380+min(energyMat.reccenterhourly.kw_energy_consumption)-10*eps(max(energyMat.reccenterhourly.kw_energy_consumption));
   pred_vals(stepNum)=x;act_vals(stepNum)=energyMat.reccenterhourly.kw_energy_consumption(stepNum+120+1);
   MSError(stepNum)=(energyMat.reccenterhourly.kw_energy_consumption(stepNum+120+1)-x);
    prevActiveSt=activeSt2;
    if length(sum(predCellsPrev))>0
    [predError,accuracy,recall,precis]=findAccuracy(activeSt2,reshape(predCellsPrev,1,cellNum*colNum));
    end
    predCellsPrev=reshape(predictedSt,cellNum,colNum);
    %Store Prediction error
    %And Accuracy for predictions in TM
    x(stepNum)=stepNum;accur(stepNum)=accuracy;anomaly_score(stepNum)=(numOfBurstingCols/numOfActiveCols);
    z(stepNum)=stepNum;pred_err(stepNum)=predError;prec(stepNum)=precis;rec(stepNum)=recall;
   %PLOTS
   if stepNum==N
        figure(1)
        plot(x,accur) 
        xlabel('Seconds')
        ylabel('Accuracy')
        %new.Prediction error
        figure(2)
        plot(z,pred_err)
        xlabel('Seconds')
        ylabel('Prediction Error')
         figure(3)
        plot(z,rec)
        xlabel('Seconds')
        ylabel('Recall')  
        figure(4)
        plot(z,prec)
        xlabel('Seconds')
        ylabel('Precision')
        figure(5)
        plot(z,anomaly_score)
        xlabel('Seconds')
        ylabel('Anomaly Score')
        
        figure(6)
        plot(z,MSError);
        xlabel('Seconds')
        ylabel('MSE error')
        figure(7)
        real_values=plot(z,act_vals);
        set(real_values,'Color','red');
        hold on
        pred_values=plot(z,pred_vals);
        set(pred_values,'Color','blue');
        xlabel('Seconds')
        ylabel('Energy Consumption value')
        hold off
        legend('Real Values in timestep t+1','Predicted values for (t+2)')
   end
end
  %End of Script
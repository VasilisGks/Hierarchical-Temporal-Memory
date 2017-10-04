function [ PredError,accuracy,recall,precis] = findAccuracy( activeSt,predictedStPrev)
%Computing accuracy , Prediction_Error for predictions
%Init
   % if length(find(predictedStPrev))>0
       TP=sum(and(activeSt,predictedStPrev));
       temp=or(activeSt,predictedStPrev);
       TN=length(find(temp==0));
       accuracy=(TP+TN)/length(activeSt);
       
       temp2=xor(activeSt,predictedStPrev);
       FP=length(find(temp2.*(activeSt==0)));
       precis=TP/(TP+FP);
       
       FN=length(find(temp2.*(activeSt==1)));
       recall=TP/(TP+FN);

       TP2=sum(and(activeSt,predictedStPrev));
       PredError=1-(TP2/length(find(activeSt)));
  %  end
end



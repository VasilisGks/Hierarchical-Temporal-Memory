classdef SDRClassifier < handle

  
properties (SetAccess= private)
  k;            % Number of timesteps into the future to predict
  alpha;        % Training rate
  predHistory;  % Keep the k latest predictions to enable online learning
  W;            % Network weights
  inputSize;
  targetSize;
end

methods
  function this= SDRClassifier(k,alpha, inputSize,targetSize)
  % inputSize: TM SDR size, targetSize: encoder buckets
    this.k= k; this.alpha= alpha;
    this.predHistory= sparse(targetSize,k);
    this.inputSize= inputSize;
    this.targetSize= targetSize;
    this.W= zeros(targetSize,inputSize);
  end
  
  function prediction= predict(this, tmsdr_t, target_t, learningOn)
    % Transform target from index to logical
    if numel(target_t)==1
      tmp= target_t;
      target_t= zeros(this.targetSize,1); target_t(tmp)= 1;
    end
    
    % Activation {sparse}
    %length(this.W)
   % length(tmsdr_t)
    act= this.W*tmsdr_t';
    % Prediction probability distribution
    prediction= exp(act) ./ sum(exp(act));
    % Adapt weights

    if learningOn
      dW= -this.alpha .* (this.predHistory(:,end)-target_t) * tmsdr_t;
      this.W= this.W + dW;
    end
    % Update prediction history
    this.predHistory= [prediction, this.predHistory(:,1:end-1)];
  end
end
end
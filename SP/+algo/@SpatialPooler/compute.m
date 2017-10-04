function active= compute(sp, x, learningOn,timestep)
% Note: rows_2dtopo= sqrt(ncolumns/aspectRatio_2dtopo) must be integer
% Arguments:
% - x [ninput]x[1] {bin}: input (encoder output)
% - sp.potentialSynapses [ncolumns]x[ninput] {uint8}: permanence forall potential synapses
% - sp.boost [ncolumns]: vector of sp.boost values for every column
% - sp.dutyCycle {struct}:
%   - .active [ncolumns]: how often the column wins (moving average)
%   - .overlap [ncolumns]: how much the column overlaps with input (moving average)
% - sp.params {struct}:
%   - .connectPermThreshold {[0-1]}: Threshold permanence, over which the synapse is connected
%   - .permMax: Permanence maximum value. If permanence is floating-point, should be 1. If uint8, 255
%   - .aspectRatio_2dtopo: In a 2-D topology, the columns are arranged in a rectangle with this aspect ratio
%   - .numActiveColumnsPerInhArea: Number of active columns in each local inhibition area
%   - .stimulusThreshold: Minimum feedforward input necessary to activate column (for noise resistance)
%   - .synPermActiveInc {[0-1]}: Amount by which positively-reinforced synapses increase
%   - .synPermInactiveDec {[0-1]}: Amount by which negatively-reinforced synapses decrease
%   - .dutyCyclePeriod: How many execution timesteps to calculate the duty cycles in
%   - .minOverlapDutyCycle

[ncolumns,ninput]= size(sp.potentialSynapses);
rows_2dtopo= sqrt(ncolumns/sp.params.aspectRatio_2dtopo);

%% 2: Compute the overlap with the current input for each column
% receptiveFields [ncolumns]x[ninput] {bin,sparse}: active feedforward connections forall columns
% overlap [ncolumns]x[1] {bin}: column overlap
receptiveFields= (sp.potentialSynapses > sp.params.connectPermThreshold*sp.params.permMax);
overlap= receptiveFields*x;
active= overlap.*sp.boost;

% Show receptive fields on input space in sequence
% figure; for i=1500:2048, imagesc(reshape(receptiveFields(i,:),16,16)'+2*reshape(x,16,16)'); title(num2str(i)); colorbar; pause(0.035); end;

%% 3: Compute the winning columns after inhibition
if sp.params.globalInhibition
  kth= sort(active,'descend');
  kth= kth(sp.params.numActiveColumnsPerInhArea);
  active= sparse(active>sp.params.stimulusThreshold & active>=kth);
else
  neighborhood= 2*floor(inhibitionRadius)+1;
  % "active" is a vector. The inhibition is based on a 2-D topology, 
  %   so "active" must be represented as a 2-D matrix here
  % TODO(optimization): replace ordfilt2 with custom from Jim's toolbox
  localkth= ordfilt2(reshape(active,rows_2dtopo,[]), ...
                     neighborhood^2+1 - sp.params.numActiveColumnsPerInhArea, ...
                     true(neighborhood),'symmetric');
  localkth= reshape(localkth,ncolumns,1);
  active= sparse(active>sp.params.stimulusThreshold & active>=localkth);
end

%% 4: Update synapse permanences and internal variables
if learningOn
  % 4.1: Update permanences
  %{
  activeIdx= find(active);  % indices of active columns
  activeSynapses= sp.potentialSynapses(activeIdx,:)>sp.params.connectPermThreshold*sp.params.permMax;  % active synapses of active columns
  % find linear indices for sp.potentialSynapses
  [i,j]= find(activeSynapses);
  actsynIdx= (activeIdx(i)-1)*ninput+1 +j;
  [i,j]= find(~activeSynapses);
  inactsynIdx= (activeIdx(i)-1)*ninput+1 +j;
  %}
  activeSynapses= false(size(sp.potentialSynapses));
  % Active synapses: belong to the potential pool (sp.potentialSynapses>0) and input there is active
  activeSynapses(active,:)= sp.potentialSynapses(active,:)>0 & repmat(x',sum(active),1);  % active synapses of active columns
  sp.potentialSynapses(activeSynapses)= min(sp.params.permMax,...  % increase active synapses
    sp.potentialSynapses(activeSynapses) + sp.params.synPermActiveInc*sp.params.permMax);
  
  inactiveSynapses= activeSynapses;
  inactiveSynapses(active,:)= sp.potentialSynapses(active,:)>0 & ~repmat(x',sum(active),1);
  sp.potentialSynapses(inactiveSynapses)= max(1,...             % decrease inactive synapses (min=1, because 0 means not in potential pool
    sp.potentialSynapses(inactiveSynapses) - sp.params.synPermInactiveDec*sp.params.permMax);
  
  %{
  % 4.2: Update sp.boosting, the permanence of low-overlap columns and inhibition radius
  % forall column neighborhoods, calculate the neighborhood's mean activeDutyCycle
  localMeanActiveDC= imboxfilt(reshape(sp.dutyCycle.active,rows_2dtopo,[]), ...  % mean filter
                               neighborhood, 'padding','symmetric');
	localMeanActiveDC= reshape(localMeanActiveDC,ncolumns,1);
  % forall columns, update {active,overlap}DutyCycle
  dcPeriod= min(sp.params.dutyCyclePeriod, timestep);
  sp.dutyCycle.active= ((dcPeriod-1)*sp.dutyCycle.active + active)./dcPeriod;
  sp.dutyCycle.overlap= ((dcPeriod-1)*sp.dutyCycle.overlap + (overlap>0))./dcPeriod;
  % forall columns, update sp.boost and increase permanences if < minOverlapDutyCycle
  sp.boost= boostFunction(sp.dutyCycle.active, localMeanActiveDC);
	sp.boost= min(sp.boost,sp.params.maxBoost);
  sp.potentialSynapses(sp.dutyCycle.overlap < sp.params.minOverlapDutyCycle,:)= ...
    sp.potentialSynapses(sp.dutyCycle.overlap < sp.params.minOverlapDutyCycle,:) + ...
    0.1*sp.params.connectPermThreshold*sp.params.permMax;
  % update inhibition radius
  inhibitionRadius= sqrt(mean(sum(receptiveFields,2)));
  %}
  %figure(2); subplot(211); histogram(sp.dutyCycle.overlap); subplot(212); imagesc(reshape(active,32,64)');
end
end

function boost= boostFunction(activity, localActivity)
boost= exp(localActivity - activity);
end

function potentialSynapses= initSynapses(sp,ninput,ncolumns,spdim,indim,mappingFraction,centerR)
  potentialSynapses= zeros(ninput,ncolumns,'uint8');
  for c= 0:ncolumns-1                                 % Foreach column...
    sppos= [mod(c,spdim(1)), floor(c/spdim(1))] +1;   % 2D index in SP-space
    [syn,center]= map_sp2in(sppos, indim,mappingFraction,centerR,sp.params.potentialDensity); % Map to input-space
    syn= pmodel(syn,center,sp.params.connectPermThreshold,sp.params.permMax);   % Generate permanences
    potentialSynapses(:,c+1)= syn;
  end
  potentialSynapses= potentialSynapses';
end

function [syn,center]= map_sp2in(sppos, indim,mappingFraction,centerR,potentialDensity)
  % Maps the position of a single minicolumn to the input space
  % syn [ninput] {bin}: vector of connected synapses
  % center: 2D index of mapped center

  center= floor([sppos(1)/mappingFraction(1), sppos(2)/mappingFraction(2)]) +1;  %assume ncolumns>=ninput
  % Generate all the points around the center within the mapped inhibition radius
  xgrid= center(1)-centerR(1) : center(1)+centerR(1);
  ygrid= center(2)-centerR(2) : center(2)+centerR(2);
  % ...and select only those that are within bounds
  [xgrid,ygrid]= meshgrid(xgrid(xgrid>0 & xgrid<=indim(1)), ygrid(ygrid>0 & ygrid<=indim(2)));
  idx= sub2ind(indim, xgrid, ygrid); idx= idx(:);
  % Randomly keep a predetermined proportion
  idx= idx( randperm(length(idx), round(potentialDensity*length(idx),0)) );

  % Connect the corresponding synapses
  syn= false(indim(1)*indim(2),1);
  syn(idx)= true;
  syn= uint8(syn);
end

function syn= pmodel(syn,center,threshold,permMax)
  % Assign random initial values to the synapses of 1 minicolumn based on a probability model
  % that rewards higher permanences to synapses towards the center
  % syn {uint8 (0|1)}

  syn(syn==1)= uint8(25*randn(length(syn(syn==1)),1,'single') + threshold*permMax);
  % TODO: higher permanences towards the center
end

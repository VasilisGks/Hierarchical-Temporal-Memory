classdef SpatialPooler < handle
properties
  params;
end

properties (Access= private)
  potentialSynapses;
  boost;
  inhibitionRadius;
  dutyCycle;
end

methods
  function sp= SpatialPooler(ninput,ncolumns,input_aspectRatio,params)
    sp.params= params;
    sp_aspectRatio= params.aspectRatio_2dtopo;
    spdim= [sp_aspectRatio*sqrt(ncolumns/sp_aspectRatio), sqrt(ncolumns/sp_aspectRatio)]; %[x,y]
    indim= [input_aspectRatio*sqrt(ninput/input_aspectRatio), sqrt(ninput/input_aspectRatio)];

    if params.globalInhibition, params.initInhibitionRadius= max(spdim); end;
    inhibitR= params.initInhibitionRadius;
    mappingFraction= spdim./indim; % how many minicolumns cover each input element in each dimension
    mappedR= round(inhibitR./mappingFraction,0);
    centerR= [2*mappedR(1)+1, 2*mappedR(2)+1];

    sp.params.numActiveColumnsPerInhArea= ceil(params.activityDensity .* ...
      min(params.initInhibitionRadius.^2, prod(spdim)));

    % Initialize the potential synapses
    sp.potentialSynapses= sp.initSynapses(ninput,ncolumns,spdim,indim,mappingFraction,centerR);

    sp.boost= ones(ncolumns,1);
    sp.dutyCycle.active= zeros(ncolumns,1); sp.dutyCycle.overlap= zeros(ncolumns,1);
  end

  active= compute(sp, x, learningOn,timestep)
end

methods (Access= private)
  potentialSynapses= initSynapses(sp, ninput,ncolumns,spdim,indim,mappingFraction,centerR);
end
end

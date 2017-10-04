function plotEncSpOut(enchist,sphist,timestep)

if timestep==1
  figure(3); clf;
else
  enchist= enchist(:,timestep-1:timestep); sphist= sphist(:,timestep-1:timestep);
  encOut= enchist(:,2); spOut= sphist(:,2);
  set(groot,'CurrentFigure',3);
  subplot(121);imagesc(reshape(encOut,25,24)');
  subplot(122);imagesc(reshape(spOut,64,32)');
end

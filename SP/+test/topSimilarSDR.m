function topsimilar= topSimilarSDR(sdr_history, k)
% Finds the top-k% most similar SDRs to the last one from the entire history.
% - sdrhistory [nSDR]x[time]
% - topsimilar {cell}:
%   - {1} [nSDR]x[k]: the selected SDRs
%   - {2} [k]: the timestep at which each one occurred

% Compute the overlap of the last SDR with each previous one
k= ceil(k/100*size(sdr_history,2)); if size(sdr_history,2)==1, k=0; end;
%k= min(k,size(sdr_history,2)-1);
overlap= sum(sdr_history(:,end) & sdr_history(:,1:end-1));
[~,topOverlapIdx]= sort(overlap,'descend');
topOverlapIdx= topOverlapIdx(1:k);
topsimilar= {sdr_history(:,topOverlapIdx), topOverlapIdx};

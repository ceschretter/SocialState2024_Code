function [bouts, lens] = detect_binarybouts_v2(binary_vector)
% DETECT_BINARYBOUTS serves to detect the start and end indices of binary events from vectors of zeros
% and ones as well as calculate the lens of these events. Importantly, it
% corrects for any spurious additional t_starts and t_ends that correspond
% with "unfinished" events at the beginning or end of the assay.

% About Version: v2 deals with ongoing events that do not
% return to 0 before assay completion.

% INPUT ARGUMENTS

%   binary_vector is a binary time-series; usually extracted from a
%   non-binary, continuous time-series with minimum and/or maximum
%   thresholding of this signal.

% OUTPUT ARGUMENTS

%   bouts is a Nx2 matrix where each row correspond to one bout, col 1 is
%   the t_start for the bout and col 2is t_end for the bout (in indices, or
%   frames for a courtship assay)
%
%   lens is a Nx1 vector where each row corresopnds to the length (in
%   frames or indices) of each bout, associated with the bouts matrix.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% detect t_start and t_end of each bout:

diff_v = diff(binary_vector); % differentiate 
v_s = find(diff_v>0); % positive slopes
v_e = find(diff_v<0); % negative slopes

% if all events are captured in their entireties, length(v_s) =
% length(v_e); however, on occasion, the assay will begin during an ongoing
% event or end during an ongoing event, leaving an "orphan" v_s or v_e at
% the beginning or end. The following corrects for this:

if length(v_s) > length(v_e) 
    bouts = [v_s(1:length(v_e)) v_e; v_s(end) length(binary_vector)]; 
end
if length(v_s) < length(v_e) 
    bouts = [v_s v_e(1:length(v_s))]; 
end
if length(v_s) == length(v_e) 
    bouts = [v_s v_e]; 
end

if sum(binary_vector) > 0
    if bouts(1,2)<bouts(1,1)
        bouts = [bouts(1:end-1,1) bouts(2:end,2)]
    end
end

% calculate bout lengths:
lens = bouts(:,2)-bouts(:,1);



end


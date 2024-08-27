function [visual_feats] = GetVisualFeatures(trk, calib, subject_ID,target_ID, samplerate, frame_range, nose, plotit)
%
% GETVISUALFEATURES calculates estimated visual features of a target fly 
% in a subject fly's visual field during free behavior assays that have been
% tracked with FlyTracker (https://github.com/kristinbranson/FlyTracker). 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXTERNAL DEPENDENCIES
%
% CircStats Statistics Toolbox should be downloaded and put in the MATLAB
% path: https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
% If the entire code library for this project is cloned, it will include
% this CircStats package.
%
%
% INPUT ARGUMENTS
%
% trk:      track.mat produced by FlyTracker, should be pre-loaded into
%           MATLAB workspace (variable type: data structure)
% calib:    calibration.mat data structure from FlyTracker, should be
%           pre-loaded in to workspace. (variable type: data structure)
% subject_ID:  Fly ID for subject fly (verify in FlyTracker visualizer; variable type: integer)
% target_ID:   Fly ID for target fly (to be considered in the visual field
%              of the subject fly; variable type: integer)
% samplerate: this is the desired sample rate; if want to produce visual features at same frame
%               that video was acquired, then samplerate = FPS of raw
%               video. But, if captured videos at different frame rate, may
%               wish to downsample to match the lowest frame rate among
%               videos in data set. Indicate this frame rate here (in Hz).
% frame_range: this is a two integer vector that designates the start and
%              end frames to be included in the analysis. If copulation
%              event occurs, may want to only include frames prior to
%              copulation start frame.  
% plotit:      true: will plot time series of all variables 
%           false: no plot will be generated 
%
%
%
% OUTPUT ARGUMENTS
%
% This script fits an ellipsoid to the target fly and produces a data structure
% containing the following time-series for the designated frame range:
%
% Time vector:          in seconds
% Angular Position:     position of target centroid with respect to subject's
%                       gaze at 0 (degrees); positive values are to the
%                       left (counter-clockwise), negative values are to the right (clockwise).
% Angular Size X:       horizontal angular diameter of target in subject's visual field (degrees)
% Angular Size Y:       vertical angular diameter of target in subject's visual field (degrees)
% Angular Area:         angular area of target in the subject's visual field (degrees^2)
% Angular Velocity:     change in angular position of target over time (degrees/second) with 
%                       respect to subject's visual field center; positive
%                       values indicate counter-clockwise movement
%                       (leftward), negative values indicate clockwise
%                       movement (rightward)              
% Angular Expansion X:  change in horizontal angular diameter of target over time (degrees/second)
% Angular Expansion Y:  change in vertical angular size of target over time (degrees/second)
% Angular Expansion Area:   change in angular area of target (degrees^2/second)
% Metric Length X:      width of subject (mm) as detected by FlyTracker
% Metric Length Y:      length of subject (mm) as detected by FlyTracker
% Metric Distance:      distance of euclidean path between subject and target centroids (mm) 
% Relative Translational Velocity:    relative velocity of target and subject = change in metric distance over time (mm/s)
%
%
% data structure also lists names of these time-series and their units 
%
% *NOTE: for variables that are derivatives, they're all calculated on a
% frame-by-frame basis. Depending on the video's frame rate, this may or
% may not be a relevant rate for this calculation. One may need to choose a
% reasonable window over which to smooth these variables. Beware that the mean and
% standard deviation will shrink over larger windows. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Downsample tracked variables if desired:

trk.data = trk.data(:,frame_range(1):frame_range(2),:);

if samplerate ~= calib.FPS
    % interpolate and downsample to desired sample rate:
    t_max = length(trk.data)./calib.FPS; % end time of video 
    original_time_vect = [1:length(trk.data)]./calib.FPS; 
    new_time_vect = [1:1:t_max*samplerate]./samplerate;
    
    [n_flies, len_trks, n_trk_vars] = size(trk.data);
    
    for i = 1:n_flies % 2 flies
        for j = 1:n_trk_vars % downsample all 35 track variables
            %new_trks.data(i,:,j) = interp1(original_time_vect, trk.data(i,:,j), new_time_vect);
            downsampled_trk(i,:,j) = interp1(original_time_vect, trk.data(i,:,j), new_time_vect);
        end
    end
    trk.data = downsampled_trk; % replace with downsampled data
    FPS = samplerate; % redefine frame rate to downsampled rate
else
    FPS = calib.FPS; % if don't downsample, the sampling rate is the same as raw video frame rate
end
    

%% Pull relevant parameters from calibration.mat

ppm= calib.PPM;
%frame_width = calib.w;
[frame_width, frame_height] = size(calib.mask);
%frame_height = calib.h;


ori1 = trk.data(subject_ID,:,3); % subject fly ori
ori2 = trk.data(target_ID,:,3); % fem ori
posx1 = trk.data(subject_ID,:,1); % x and y coordinates of subject fly centroid
posy1 = trk.data(subject_ID,:,2);
posx2 = trk.data(target_ID,:,1); % x and y coordinates of target fly centroid
posy2 = trk.data(target_ID,:,2);
major_axis_vect = trk.data(target_ID,:,4)./ppm; % major ellipse axis
minor_axis_vect = trk.data(target_ID,:,5)./ppm; % minor ellipse axis
   
% transpose coordinates so that (0,0) is equiavlent to bottom left corner of frame, and therefore all coordinates are positive values:

center_m = [posx1' frame_height-posy1'];
center_f = [posx2' frame_height-posy2']; 

display('Calculating visual parameters for all frames ');

%% Translate + rotate fly centroids so subject is at the origin:

trans_ctr_subject = [center_m(:,1)-center_m(:,1) center_m(:,2)-center_m(:,2)]; % should be [x y] = [0 0] for all frames
trans_ctr_target = [center_f(:,1)-center_m(:,1) center_f(:,2)-center_m(:,2)];  % should be same distance and angle from subject fly for every frame, but origin = subject fly 

% rotate fly centers so that subject fly is facing 0 radians = 0 deg:
theta = -1*ori1' + deg2rad(0);

for i = 1:length(ori1') 
    clear rot_temp;
    rot_temp = [cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))]; % in radians!
    rotated_f(i,:) = rot_temp*(trans_ctr_target(i,:)');
    new_ori_f(i,:) = theta(i) + ori2(i);
    theta_target(i,:) = cart2pol(rotated_f(i,1),rotated_f(i,2));
end

theta_f2m = deg2rad(theta_target)';
ctr_target = rotated_f; % cartesian coordinates of target fly relative to subject fly in new basis

%% Fit ellipse to target fly and generate points around perimeter:

[xe,ye] = get_ellipse(trk, ctr_target, new_ori_f, target_ID);


%% If want to calculate everything with respect to subject fly's "nose" (anterior-most point), rather than centroid:

[m,n] = size(xe); % m = number of frames, n = number of points ellipse fit perimeter

if nose
    for i = 1:m
        xe(i,:) = xe(i,:)-(trk.data(subject_ID,i,4))/2;
    end
end



%% Calculate linear fit slope (m) between subject fly and target fly in this new basis:
L = []; % intitializing
x = ctr_target(:,1); 
y = ctr_target(:,2); 

for i = 1:length(x) %across every frame, calculate m_norm and offset value for translation (point b)
    m(i) = (0 - y(i))/(0 - x(i));
    m_norm(i) = -(1/m(i));  % perpendicular slope to line from target fly to subject fly is opposite reciprocol
    % determine y-intercept of fxn that passes through target fly centroid with
    % slope m_norm
    b(i) = y(i)-m_norm(i)*(x(i)); % b parameter for function that is normal to line between target fly and subject fly
    % y = m_norm*x+b
    xintercepts(i) = (0-b(i))/m_norm(i); %intercept (offset for translation) is the pivot point for rotation (theta_rot)
end

%% Translate and rotate target ellipse perimeter points and calculate output variables across all frames:

for i = 1:length(ctr_target) % across all frames!

if isnan(ctr_target(i,:)) == 0    
    clear L_t
    % translation of ellipse perimeter points so that pivot point (b) is at origin for rotation:
    offset = xintercepts(i); %tan(2*pi-theta_f2m(i)).*ctr_target(i,2)+ctr_target(i,1); 
    trans_ellipse = [xe(i,:)'-offset  ye(i,:)']; % translate so rotation center is at origin
    
    %rotation about xintercept (new origin following translation) of line normal to line between subject fly/target fly
    theta_rot = -(theta_f2m(i)-1.5*pi); %-atan(trans_centroid(2)/trans_centroid(1)); % rotation angle to put target fly ellipse on x-axis

    %logical arguments to determine if rotation is positive or negative theta_rot:
    
    if x>0 & y>0
        theta_rot = theta_rot; % rotation angle to put target fly ellipse on x-axis
    end
    if x<0 & y<0
        theta_rot = theta_rot;
    end
    if x>0 & y<0
        theta_rot = -theta_rot;
    end
    if x<0 & y>0
        theta_rot = -theta_rot;
    end
    
    rot_mat = [cos(theta_rot) -sin(theta_rot); sin(theta_rot) cos(theta_rot)]; % in radians!
    
    for j = 1:length(trans_ellipse) %across all points in ellipse 
        rotated_ellipse(j,:) = rot_mat*(trans_ellipse(j,:)');
    end
   
%% Do you want visual features with respect to subject fly's anterior-most point ("nose")?
%   If so, subtract half of subject fly body length from target fly x
%   coordinates:

  




%% Calculate Output Variables across all frames from rotated and translated target ellipse:

    max_x = max(rotated_ellipse(:,1));
    min_x = min(rotated_ellipse(:,1));
    max_ind = find(rotated_ellipse(:,1) == max_x); %max(xp)
    min_ind = find(rotated_ellipse(:,1) == min_x); %min(xp) as viewed in subject fly VF
    P1 = rotated_ellipse(max_ind,:); P2 = rotated_ellipse(min_ind,:);
    L_t = pdist([P1; P2]); %true length of target fly object in subject fly VF
    % in original subject fly basis, re-locate these same indices (points in
    % ellipse)
    max_coords = [xe(i,max_ind) ye(i,max_ind)];
    min_coords = [xe(i,min_ind) ye(i,min_ind)];
    midpt_coords = (min_coords + max_coords).'/2; % real centroid in subject fly VF
   
    % define x and y vals for logical arguments (assignment of radian vals)
    x =  midpt_coords(1); y =  midpt_coords(2);
    x_max = max_coords(1); y_max = max_coords(2);
    x_min = min_coords(1); y_min = min_coords(2);
   
    %logical arguments:

    if x>0 & y>0
        ang_pos_ctr = atan(y/x);
        ang_pos_max = atan(y_max/x_max);
        ang_pos_min = atan(y_min/x_min);
    end
    if x<0 & y>0
        ang_pos_ctr = pi+atan(y/x);
        ang_pos_max = pi+atan(y_max/x_max);
        ang_pos_min = pi+atan(y_min/x_min);
    end
    if x<0 & y<0
        ang_pos_ctr = pi+atan(y/x);
        ang_pos_max = pi+atan(y_max/x_max);
        ang_pos_min = pi+atan(y_min/x_min);
    end
    if x>0 & y<0
        ang_pos_ctr = 2*pi+atan(y/x);
        ang_pos_max = 2*pi+atan(y_max/x_max);
        ang_pos_min = 2*pi+atan(y_min/x_min);
    end
    
  
   % Calculate position and size parameters for each frame:
   ang_centroid(i,1) = wrapToPi(ang_pos_ctr); % true angular position of target fly object in subject fly VF (rad)  
   lin_dist(i,1) = pdist([x y ; 0 0]); % euclidean distance between origin (subject fly) and target fly object as viewed in subject fly VF (pixels)
   ang_size_x(i,1) = (atan((0.5*L_t)/lin_dist(i))); % half-size in radians
   ang_size_y(i,1) = (atan((0.5*ppm)/lin_dist(i))); ; % assume constant target fly posture in z dimension of ~ 1 mm (half-size in radians)
   ang_size_x(i,1) = 2*rad2deg(ang_size_x(i,1)); % full size in degrees
   ang_size_y(i,1) = 2*rad2deg(ang_size_y(i,1)); % full size in degrees
   ang_area(i,1) = (ang_size_x(i,1)/2)*(ang_size_y(i,1)/2)*pi; %(ellipse area = a*b*pi; in deg^2)
   L = [L; L_t]; % "real" length of animal in pixels 
   
   
else
    
   ang_centroid(i,1) = NaN; % true angular position of target fly object in subject fly VF (rad)
   lin_dist(i,1) = NaN; % euclidean distance between origin (subject fly) and target fly object as viewed in subject fly VF (pixels)
   ang_size_x(i,1) = NaN; % half-size in radians
   ang_size_y(i,1) = NaN; % assume constant target fly posture in z dimension of ~ 1 mm (half-size in radians)
   ang_area(i,1) = NaN; %(ellipse area = a*b*pi; in deg^2)
   L = [L; NaN]; % in pixels for now
end
    
end

debugpoint = 1; 


% Linearly interpolate any NaNs:
ang_size_x = fillmissing(ang_size_x,'linear');
ang_size_y = fillmissing(ang_size_y,'linear');
ang_area = fillmissing(ang_area,'linear');

% calculate derivatives of non-circular variables (frame-by-frame):
diff_ang_size_x = diff(ang_size_x).*FPS; 
diff_ang_size_y = diff(ang_size_y).*FPS;
diff_ang_area = diff(ang_area).*FPS;
lin_vel = diff(lin_dist).*FPS;

% calculate angular velocity using circ stats toolbox:

for i = 1:length(ang_centroid)-1
    if sum(isnan(ang_centroid([i i+1]))) == 0
    ang_vel(i,1) = circ_dist2(ang_centroid(i+1), ang_centroid(i));
    else
    ang_vel(i,1) = NaN;   
    end
end
ang_vel = rad2deg(ang_vel)*FPS; % convert to deg/sec
ang_vel = smooth(ang_vel,3); % this vector can be a bit noisy without any smoothing

% pad all with one NaN:
diff_ang_size_x = [NaN; diff_ang_size_x];
diff_ang_size_y = [NaN; diff_ang_size_y];
diff_ang_area = [NaN; diff_ang_area];
lin_vel = [NaN; lin_vel];
ang_vel = [NaN; ang_vel];

% interpolate for later use (won't be using very first frame anyway):
diff_ang_size_x = fillmissing(diff_ang_size_x,'linear');
diff_ang_size_y = fillmissing(diff_ang_size_y,'linear');
diff_ang_area = fillmissing(diff_ang_area,'linear');
lin_vel = fillmissing(lin_vel,'linear');
ang_vel = fillmissing(ang_vel,'linear');

%% Get rid of spurious misdetections of angular position (due to transient orientation flips of subject fly):
% detect with high angular velocities of target fly

inds = find(abs(ang_vel)>500);
ang_vel(inds)=NaN;
ang_vel = fillmissing(ang_vel,'linear');
ang_centroid(inds)=NaN;
ang_centroid = fillmissing(ang_centroid,'linear');

%% Output Data Structure:

% Time series of calculated visual parameters:
visual_feats.data(:,1) = (1:length(ang_centroid))./FPS; % time vector (sec)
visual_feats.data(:,2) = rad2deg(ang_centroid); % true centroid of target fly in subject fly VF (deg), in x
visual_feats.data(:,3) = ang_size_x; % angular size in horizontal (deg)
visual_feats.data(:,4) = ang_size_y; % angular size in vertical (deg)
visual_feats.data(:,5) = ang_area; % angular area of ellipse fit to major and minor axes of ang_size_x and ang_size_y, respectively (deg^2)
visual_feats.data(:,6) = ang_vel; % angular velocity of true target fly centroid in subject fly VF; change in ang_centroid over time (deg/sec)
visual_feats.data(:,7) = diff_ang_size_x; %change in horizontal angular size over time (deg/sec)
visual_feats.data(:,8) = diff_ang_size_y; %change in vertical angular size over time (deg/sec)
visual_feats.data(:,9) = diff_ang_area; %change in angular area over time (deg^2/s)
visual_feats.data(:,10) = L./ppm; % true horizontal length (in x) of target fly in subject fly VF (mm) 
visual_feats.data(:,11) = 1; % approximate vertical length (in y) of target fly in subject fly VF (mm) 
visual_feats.data(:,12) = lin_dist./ppm; % distance of euclidean path between subject fly and target fly centroid in subject fly VF (mm)
visual_feats.data(:,13) = lin_vel./ppm; %relative velocity of target fly and subject fly = change in lin_dist over time


% Name for each time series:
visual_feats.names(1) = {'Time'};
visual_feats.names(2) = {'Angular Position'};
visual_feats.names(3) = {'Angular Size X'};
visual_feats.names(4) = {'Angular Size Y'};
visual_feats.names(5) = {'Angular Area'};
visual_feats.names(6) = {'Angular Velocity'};
visual_feats.names(7) = {'Diff Angular Size X'};
visual_feats.names(8) = {'Diff Angular Size Y'};
visual_feats.names(9) = {'Diff Angular Area'};
visual_feats.names(10) = {'Metric Length X'};
visual_feats.names(11) = {'Metric Length Y'};
visual_feats.names(12) = {'Metric Distance'};
visual_feats.names(13) = {'Metric Relative Velocity'};


% Units of each time series:
visual_feats.units(1) = {'sec'};
visual_feats.units(2) = {'degrees'};
visual_feats.units(3) = {'degrees'};
visual_feats.units(4) = {'degrees'};
visual_feats.units(5) = {'degrees^2'};
visual_feats.units(6) = {'degrees/second'};
visual_feats.units(7) = {'degrees/second'};
visual_feats.units(8) = {'degrees/second'};
visual_feats.units(9) = {'degrees^2/second'};
visual_feats.units(10) = {'millimeters'};
visual_feats.units(11) = {'millimeters'};
visual_feats.units(12) = {'millimeters'};
visual_feats.units(13) = {'millimeters/second'};


%save([prefix '-visual_feats'],'visual_feats')
display(['Visual features calculated!']) % and saved in ' num2str([cd]) '!']);

%% Optional: If want to plot some time series and histograms, uncomment below:

if plotit
    figure
    vectors_to_plot = [2 3 4 5 6 12];
    for i = 1:6
        subplot(6,1,i)
        plot(visual_feats.data(:,1), visual_feats.data(:,vectors_to_plot(i)),'k', 'LineWidth',2)
        %ylabel([visual_feats.names(vectors_to_plot(i)) ' (' visual_feats.units(vectors_to_plot(i)) ')'],'FontSize',12)
        xlabel('Time (s)','FontSize',12)
        box off
        xlim([0 300])
    end

    subplot(6,1,1); ylabel(['Position (' char(176) ')']); title('Basic Visual Features of Target Fly')
    subplot(6,1,2); ylabel(['Diameter in X (' char(176) ')'])
    subplot(6,1,3); ylabel(['Diameter in Y (' char(176) ')'])
    subplot(6,1,4); ylabel(['Area (deg^2)'])
    subplot(6,1,5); ylabel(['Velocity (' char(176) '/s)']); ylim([-500 500])
    subplot(6,1,6); ylabel(['Distance (mm)'])

    figure
    for i = 1:6
        subplot(6,1,i)
        histogram(visual_feats.data(:,vectors_to_plot(i)), 500,'FaceColor','k')
        ylabel('Frames')
        box off
    end
    subplot(6,1,1); xlabel(['Position (' char(176) ')']); title('Basic Visual Features of Target Fly'); xlim([-60 60])
    subplot(6,1,2); xlabel(['Diameter in X (' char(176) ')']); xlim([0 120])
    subplot(6,1,3); xlabel(['Diameter in Y (' char(176) ')']); xlim([0 45])
    subplot(6,1,4); xlabel(['Area (deg^2)']); xlim([0 1200])
    subplot(6,1,5); xlabel(['Velocity (' char(176) '/s)']); xlim([-500 500])
    subplot(6,1,6); xlabel(['Distance (mm)']); xlim([0 30])


end



end




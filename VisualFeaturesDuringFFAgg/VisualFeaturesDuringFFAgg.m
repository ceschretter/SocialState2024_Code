%% This is for examining the visual features flies experience during aggression
% What visual stimuli do the female flies experience during an aggressive interaction?

%Add paths to Male Courtship Feature Analysis
addpath '/GitHub/SocialState2024_private/Code/VisualFeaturesDuringCourtship'
addpath '/GitHub/SocialState2024_private/Code/VisualFeaturesDuringCourtship/CircStat2012a'
addpath '/GitHub/SocialState2024_private/Code/VisualFeaturesDuringCourtship/Colormaps/Colormaps (5)/Colormaps'
addpath '/GitHub/SocialState2024_private/Code/VisualFeaturesDuringCourtship/Functions'

%Add paths to data
addpath '/GitHub/SocialState2024_private/Data/FFAggression'

% Setup the calibration file, experimental directories, trajectories, JAABA files, and perframe features for experiments

calib = load('/GitHub/SocialState2024_private/Data/FFAggression/aIPg_calib.mat');
tracks = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/registered_trx.mat');
expdir = '/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR';

behaviors = load('/groups/rubin/home/schretterc/Documents/VisionProject_Code/SampleExpts/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/scores_FeAggfRGBub.mat');

angle = load(fullfile(expdir,'/','perframe','angleonclosestfly.mat'),'data');
closestID = load(fullfile(expdir,'/','perframe','closestfly_center.mat'),'data');

%% Pull parameters from -calibration.mat
FPS = calib.calib.FPS;
ppm= calib.calib.PPM;
frame_width = calib.calib.w;
frame_height = calib.calib.h;

%% Determine which frames to consider across video based on input aggression classifier:
nflies = length(size(tracks.trx.id,1)); %number of flies with trajectories 

for f = 1:nflies
    % define subject and get aggression data for subject
    subject_ind = f;
    aggflyi = behaviors.allScores.postprocessed{subject_ind};
    substart = tracks.trx(subject_ind).firstframe;
    aggflyi = aggflyi(:,substart:end);

    %compute values for the subject...
    ori1 = (tracks.trx(subject_ind).theta(:)).'; % ori of fly 1
    posx1 = (tracks.trx(subject_ind).x(:)).'; % x and y coordinates of fly1 center
    posy1 = (tracks.trx(subject_ind).y(:)).';
    posa1 = (tracks.trx(subject_ind).a(:)).';
    posn1 = posy1 - posa1;
    nose_m = [posx1' frame_height-posn1']; % nose of fly1 (furthest point on the ellipse)
    center_m = [posx1' frame_height-posy1'];

    dur = tracks.trx(subject_ind).nframes; 
    display('Calculating Visual Parameters for all frames');

    %compute values for target but only for closestID while aggressive
    ori2 = nan(size(aggflyi));
    posx2 = nan(size(aggflyi));
    posy2 = nan(size(aggflyi));
    major_axis_vect = nan(size(aggflyi));
    minor_axis_vect = nan(size(aggflyi));

    %define values used for fiting the female ellipse here
    theta = -1*ori1' + deg2rad(0); 
    phi = 0:0.1:2*pi;
    xe = zeros(length(theta), length(phi)); ye = zeros(length(theta), length(phi));
    ctr_f = [];
    ctr_coords = [];
 
    % find the closest fly perframe for that subject
    targetvalues = cell2mat(closestID.data(subject_ind));

    tic %starting timer on calculating for fly i

    for k = 1:numel(aggflyi)
        if aggflyi(k) == 1 
            target_ind = targetvalues(k);
            %organizing based on frame number 
            if substart == 1
                t = k + tracks.trx(target_ind).off ;
            else
                t = (substart + (k-1)) + tracks.trx(target_ind).off;
            end

            ori2(k) = tracks.trx(target_ind).theta(t); % ori of fly 2 
            posx2(k) = tracks.trx(target_ind).x(t); % x and y coordinates of fly2 center
            posy2(k) = tracks.trx(target_ind).y(t);
            major_axis_vect(k) = (tracks.trx(target_ind).a(t)*2)./ppm; % major ellipse axis (changed based on Alice's notes that in their code this is slightly different)
            minor_axis_vect(k) = (tracks.trx(target_ind).b(t)*2)./ppm; % minor ellipse axis 
            
            center_f = [posx2' frame_height-posy2'];
            
            % Translate + rotate fly centroids so female 1 is at origin:
            trans_ctr_m = [nose_m(k,1)-nose_m(k,1) nose_m(k,2)-nose_m(k,2)]; % should be [x y] = [0 0] for all frames
            trans_ctr_f = [center_f(k,1)-nose_m(k,1) center_f(k,2)-nose_m(k,2)];  % should be same distance and angle from male for every frame, but origin = male nose 

            % rotate fly centers so that male is facing 0 radians = 0 deg:
            theta = -1*ori1(k) + deg2rad(0);

            for i = 1:length(k) 
                clear rot_temp;
                rot_temp = [cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))]; % in radians!
                rotated_f(i,:) = rot_temp*(trans_ctr_f(i,:)');
                new_ori_f(i,:) = theta(i) + ori2(i);
                theta_target(i,:) = cart2pol(rotated_f(i,1),rotated_f(i,2));
            end

            theta_f2m = deg2rad(theta_target)';
            ctr_f_k = rotated_f; %cartesian coordinates of female relative to male in new basis

            
            %define fit the ellipse to female and generate points around
            %in frame t
            ctr_coords_k = ctr_f_k;
            theta = new_ori_f;
            phi = 0:0.1:2*pi;
            costheta = cos(theta);
            sintheta = sin(theta);


            x = ctr_coords_k(:,1); % fly x position
            y = ctr_coords_k(:,2); % fly y position
            major_axis_len = (tracks.trx(target_ind).a(t)*2); % fly major axis length 
            minor_axis_len = (tracks.trx(target_ind).b(t)*2); % fly minor axis length 
            a = major_axis_len./2;
            b = minor_axis_len./2;
            
            for i = length(t)
                X1 = a(1,i)*cos(phi);
                Y1 = b(1,i)*sin(phi);
                rot_mat = [costheta(i) -sintheta(i); sintheta(i) costheta(i)];
                    for j = 1:length(phi)
                        rot_ellipse(j,:) = rot_mat*([X1(1,j) Y1(1,j)]');
                    end
            end 
            xe(k,:) = rot_ellipse(:,1)' + x(i);
            ye(k,:) = rot_ellipse(:,2)' + y(i);
            ctr_f(k,:) = ctr_f_k;
            ctr_coords(k,:) = ctr_coords_k;
            
        else 
            target_ind = targetvalues(k);
            t = tracks.trx(subject_ind).firstframe + tracks.trx(target_ind).off;

            ori2(k) = NaN; % ori of fly 1
            posx2(k) = NaN; % x and y coordinates of fly2 center
            posy2(k) = NaN;
            major_axis_vect(k) = NaN; % major ellipse axis
            minor_axis_vect(k) = NaN;
            xe(k,:) = NaN;
            ye(k,:) = NaN;
            ctr_f(k,:) = NaN(1,2);
            ctr_coords(k,:) = NaN(1,2);

        end
        clear target_ind; %clearing to compute for next frame
    end

    rotated_f = ctr_f;

    toc % timing since started calculation

    display('Female elliptical fits complete!');
    %% Calculate linear fit slope (m) between female subject and female target in this new basis:
    %clear before start next set
    clear a b angle costheta ctr_coords ctr_coords_k ctr_f_k major_axis_len minor_axis_len phi rot_ellipse rot_mat sintheta t x y X1 Y1;
    
    [m,n] = size(xe); % m = number of frames, n = number of points ellipse fit perimeter
    L = []; % intitializing
    % Initialize (make empty values for later)
    ang_centroid = [];
    ang_size_x = [];
    ang_size_y = [];
    ang_area = [];
    lin_dist = [];
    lin_vel = [];
    x = ctr_f(:,1); 
    y = ctr_f(:,2); 

    for i = 1:length(x) %across every frame, calculate m_norm and offset value for translation (point b)
        m(i) = (0 - y(i))/(0 - x(i));
        m_norm(i) = -(1/m(i));  % perpendicular slope to line from female to male is opposite reciprocol
        % determine y-intercept of fxn that passes through female centroid with
        % slope m_norm
        b(i) = y(i)-m_norm(i)*(x(i)); % b parameter for fxn that is normal to line between female and male
        % y = m_norm*x+b
        xintercepts(i) = (0-b(i))/m_norm(i); %intercept (offset for translation) IS the pivot point for theta_rot
    end

    %% Translate + Rotate Female Ellipse Perimeter Points

    for i = 1:length(ctr_f) %across all frames!
        p = sum(xe(i,:));
        q = ctr_f(i,1);
        r = ctr_f(i,2);

        if all((~isnan(q)) & (~isnan(r)) & (p ~= 0) & (~isnan(p))) %added in another qualifier if row was zero   
            clear L_t
        % translation of ellipse perimeter points so that pivot point (b) is at origin for rotation:
    
            offset = xintercepts(i); %tan(2*pi-theta_f2m(i)).*ctr_f(i,2)+ctr_f(i,1); 
            trans_ellipse = [xe(i,:)'-offset  ye(i,:)']; % translate so rotation center is at origin
 
    % rotation about xintercept (new origin following translation) of line normal to line between male/female

            theta_rot = -(theta_f2m-1.5*pi); %-atan(trans_centroid(2)/trans_centroid(1)); % rotation angle to put female ellipse on x-axis

        %logical arguments determine if rotation is positive or negative
        %theta_rot:
    
            if x>0 & y>0
                theta_rot = theta_rot; % rotation angle to put female ellipse on x-axis
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
    


    %% Calculate Output Variables from Rotated + Translated Female Ellipse:

            max_x = max(rotated_ellipse(:,1));
            min_x = min(rotated_ellipse(:,1));
            max_ind = find(rotated_ellipse(:,1) == max_x); %max(xp)
            min_ind = find(rotated_ellipse(:,1) == min_x); %min(xp) as viewed in male VF
            P1 = rotated_ellipse(max_ind,:); P2 = rotated_ellipse(min_ind,:);
            L_t = pdist([P1; P2]); %true length of female object in male VF
            % in original male basis, re-locate these same indices (points in
            % ellipse)
            max_coords = [xe(i,max_ind) ye(i,max_ind)];
            min_coords = [xe(i,min_ind) ye(i,min_ind)];
            midpt_coords = (min_coords + max_coords).'/2; % real centroid in male VF
   
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
            ang_centroid(i,1) = wrapToPi(ang_pos_ctr); % true angular position of female object in male VF (rad)  
            lin_dist(i,1) = pdist([x y ; 0 0]); % euclidean distance between origin (male) and female object as viewed in male VF (pixels)
            ang_size_x(i,1) = (atan((0.5*L_t)/lin_dist(i))); % half-size in radians
            ang_size_y(i,1) = (atan((0.5*ppm)/lin_dist(i))); ; % assume constant female posture in z dimension of ~ 1 mm (half-size in radians)
            ang_size_x(i,1) = 2*rad2deg(ang_size_x(i,1)); % full size in degrees
            ang_size_y(i,1) = 2*rad2deg(ang_size_y(i,1)); % full size in degrees
            ang_area(i,1) = (ang_size_x(i,1)/2)*(ang_size_y(i,1)/2)*pi; %(ellipse area = a*b*pi; in deg^2)
            L = [L; L_t]; % in pixels for now
   
   
        else
    
            ang_centroid(i,1) = NaN; % true angular position of female object in male VF (rad)
            lin_dist(i,1) = NaN; % euclidean distance between origin (male) and female object as viewed in male VF (pixels)
            ang_size_x(i,1) = NaN; % half-size in radians
            ang_size_y(i,1) = NaN; % assume constant female posture in z dimension of ~ 1 mm (half-size in radians)
            ang_size_x(i,1) = NaN; % full size in degrees
            ang_size_y(i,1) = NaN; % full size in degrees
            ang_area(i,1) = NaN; %(ellipse area = a*b*pi; in deg^2)
            L = [L; NaN]; % in pixels for now
        end

            clear p v

    end

    % calculate derivatives of non-circular variables (frame-by-frame)
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
    
    % Do padding but no interpolation
    % pad all with one NaN:
    diff_ang_size_x = [NaN; diff_ang_size_x];
    diff_ang_size_y = [NaN; diff_ang_size_y];
    diff_ang_area = [NaN; diff_ang_area];
    lin_vel = [NaN; lin_vel];
    ang_vel = [NaN; ang_vel];

    % Get rid of spurious misdetections of angular position (due to transient orientation flips of subject fly):
    % detect with high angular velocities of target fly
    inds = find(abs(ang_vel)>500);
    ang_vel(inds)=NaN;
%     ang_vel = fillmissing(ang_vel,'linear');
    ang_centroid(inds)=NaN;
%     ang_centroid = fillmissing(ang_centroid,'linear');

%%

    % Output Data Structure (in desired units!):
    % Time-Series of Calculated Visual Parameters:
        visual_params.data(:,1) = (1:length(ang_centroid))./FPS; % time vector (sec)
        visual_params.data(:,2) = rad2deg(ang_centroid); % true centroid of female in male VF (deg), in x
        visual_params.data(:,3) = ang_size_x; % angular size in horizontal (deg)
        visual_params.data(:,4) = ang_size_y; % angular size in vertical (deg)
        visual_params.data(:,5) = ang_area; % angular area of ellipse fit to major and minor axes of ang_size_x and ang_size_y, respectively (deg^2)
        visual_params.data(:,6) = ang_vel; % angular velocity of true female centroid in male VF; change in ang_centroid over time (deg/sec)
        visual_params.data(:,7) = diff_ang_size_x; %change in horizontal angular size over time (deg/sec)
        visual_params.data(:,8) = diff_ang_size_y; %change in vertical angular size over time (deg/sec)
        visual_params.data(:,9) = diff_ang_area; %change in angular area over time (deg^2/s)
        visual_params.data(:,10) = L./ppm; % true horizontal length (in x) of female in male VF (mm) 
        visual_params.data(:,11) = 1; % approximate vertical length (in y) of female in male VF (mm) 
        visual_params.data(:,12) = lin_dist./ppm; % distance of euclidean path between male and female centroid in male VF (mm)
        visual_params.data(:,13) = lin_vel./ppm; %relative velocity of female and male = change in lin_dist over time
        visual_params.data(:,14) = subject_ind;


        % Name for each time-series, subject_ind, (also maybe want target_ind used..) :
        visual_params.names(1) = {'Time Vector'};
        visual_params.names(2) = {'Angular Position'};
        visual_params.names(3) = {'Angular Size X'};
        visual_params.names(4) = {'Angular Size Y'};
        visual_params.names(5) = {'Angular Area'};
        visual_params.names(6) = {'Angular Velocity'};
        visual_params.names(7) = {'Diff Angular Size X'};
        visual_params.names(8) = {'Diff Angular Size Y'};
        visual_params.names(9) = {'Diff Angular Area'};
        visual_params.names(10) = {'Metric Length X'};
        visual_params.names(11) = {'Metric Length Y'};
        visual_params.names(12) = {'Metric Distance'};
        visual_params.names(13) = {'Metric Relative Velocity'};
        visual_params.names(14) = {'Fly examining'};



        % Units of each time series:
        visual_params.units(1) = {'seconds'};
        visual_params.units(2) = {'degrees'};
        visual_params.units(3) = {'degrees'};
        visual_params.units(4) = {'degrees'};
        visual_params.units(5) = {'degrees^2'};
        visual_params.units(6) = {'degrees/second'};
        visual_params.units(7) = {'degrees/second'};
        visual_params.units(8) = {'degrees/second'};
        visual_params.units(9) = {'degrees^2/second'};
        visual_params.units(10) = {'millimeters'};
        visual_params.units(11) = {'millimeters'};
        visual_params.units(12) = {'millimeters'};
        visual_params.units(13) = {'millimeters/second'};
        visual_params.units(14) = {'fly number in tracks'};

        %save visual params for this fly (subject_ind)
        fullname = ['flysubject' num2str(subject_ind) '-visual_params.mat'];
        fullpath = fullfile(expdir, fullname);
        save(fullpath, 'visual_params');
        display(['Visual parameters calculated and saved to:' fullpath]);
        % clear varibles going to use for next fly...
        clear all

        %reload the key varibles needed
        calib = load('/GitHub/SocialState2024_private/Data/FFAggression/aIPg_calib.mat');
        tracks = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/registered_trx.mat');
        expdir = '/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR';
        behaviors = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/scores_FeAggfRGBub.mat');
        angle = load(fullfile(expdir,'/','perframe','angleonclosestfly.mat'),'data');
        closestID = load(fullfile(expdir,'/','perframe','closestfly_center.mat'),'data');
        FPS = calib.calib.FPS;
        ppm= calib.calib.PPM;
        frame_width = calib.calib.w;
        frame_height = calib.calib.h;
        nflies = length(size(tracks.trx.id,1));

end

%% Plot all visual params together
addpath '/GitHub/SocialState2024_private/Code/VisualFeaturesDuringCourtship/Functions';

% Setup the calibration file, trajectories, and experimental dirs being used 
calib = load('/GitHub/SocialState2024_private/Data/FFAggression/aIPg_calib.mat');
tracks = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/registered_trx.mat');
expdir = '/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR';

FPS = calib.calib.FPS;
ppm= calib.calib.PPM;
frame_width = calib.calib.w;
frame_height = calib.calib.h;

% Load in all of the visual params you just calculated below here:
fly1 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject1-visual_params.mat');
fly2 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject2-visual_params.mat');
fly3 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject3-visual_params.mat');
fly4 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject4-visual_params.mat');
fly5 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject5-visual_params.mat');
fly6 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject6-visual_params.mat');
fly7 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject7-visual_params.mat');
fly8 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject8-visual_params.mat');
fly9 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject9-visual_params.mat');
fly10 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject10-visual_params.mat');
fly11 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject11-visual_params.mat');
fly12 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject12-visual_params.mat');
fly13 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject13-visual_params.mat');
fly14 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject14-visual_params.mat');
fly15 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject15-visual_params.mat');
fly16 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject16-visual_params.mat');
fly17 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject17-visual_params.mat');
fly18 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject18-visual_params.mat');
fly19 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject19-visual_params.mat');
fly20 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject20-visual_params.mat');
%Fly21 was not used because it was all NaNs
fly22 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject22-visual_params.mat');
fly23 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject23-visual_params.mat');
fly24 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject24-visual_params.mat');
fly25 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject25-visual_params.mat');
fly26 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject26-visual_params.mat');
fly27 = load('//GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject27-visual_params.mat');
fly28 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject28-visual_params.mat');
fly29 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject29-visual_params.mat');
fly30 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject30-visual_params.mat');
fly31 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject31-visual_params.mat');
fly32 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject32-visual_params.mat');
fly33 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject33-visual_params.mat');
fly34 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject34-visual_params.mat');
fly35 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject35-visual_params.mat');
fly36 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject36-visual_params.mat');
fly37 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject37-visual_params.mat');
fly38 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject38-visual_params.mat');
fly39 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject39-visual_params.mat');
fly40 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject40-visual_params.mat');
fly41 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject41-visual_params.mat');
fly42 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject42-visual_params.mat');
fly43 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject43-visual_params.mat');
fly44 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject44-visual_params.mat');
fly45 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject45-visual_params.mat');
fly46 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject46-visual_params.mat');
fly47 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject47-visual_params.mat');
fly48 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject48-visual_params.mat');
fly49 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject49-visual_params.mat');
fly50 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject50-visual_params.mat');
fly51 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject51-visual_params.mat');
fly52 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject52-visual_params.mat');
fly53 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject53-visual_params.mat');
fly54 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject54-visual_params.mat');
fly55 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject55-visual_params.mat');
fly56 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject56-visual_params.mat');
fly57 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject57-visual_params.mat');
fly58 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject58-visual_params.mat');
fly59 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject59-visual_params.mat');
fly60 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject60-visual_params.mat');
fly61 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject61-visual_params.mat');
fly62 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject62-visual_params.mat');
fly63 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject63-visual_params.mat');
fly64 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject64-visual_params.mat');
fly65 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject65-visual_params.mat');
fly66 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject66-visual_params.mat');
fly67 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject67-visual_params.mat');
fly68 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject68-visual_params.mat');
fly69 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject69-visual_params.mat');
fly70 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject70-visual_params.mat');
fly71 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject71-visual_params.mat');
fly72 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject72-visual_params.mat');
fly73 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject73-visual_params.mat');
fly74 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject74-visual_params.mat');
fly75 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject75-visual_params.mat');
fly76 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject76-visual_params.mat');
fly77 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject77-visual_params.mat');
fly78 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject78-visual_params.mat');
fly79 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject79-visual_params.mat');
fly80 = load('/GitHub/SocialState2024_private/Data/FFAggression/20231107T081717_rig1_flyBubble1_SSempty_72C11LexAivk5LexAopTrpAVIE260UASGtACR_20230419_AggAVLPGtACR/flysubject80-visual_params.mat');


alldata = [fly1.visual_params.data;fly2.visual_params.data;fly3.visual_params.data;fly4.visual_params.data;fly5.visual_params.data;...
    fly6.visual_params.data;fly7.visual_params.data;fly8.visual_params.data;fly9.visual_params.data;fly10.visual_params.data;...
    fly11.visual_params.data;fly12.visual_params.data;fly13.visual_params.data;fly14.visual_params.data;fly15.visual_params.data;...
    fly16.visual_params.data;fly17.visual_params.data;fly18.visual_params.data;fly19.visual_params.data;fly20.visual_params.data;...
    fly22.visual_params.data;fly23.visual_params.data;fly24.visual_params.data;fly25.visual_params.data;...
    fly26.visual_params.data;fly27.visual_params.data;fly28.visual_params.data;fly29.visual_params.data;fly30.visual_params.data;...
    fly31.visual_params.data;fly32.visual_params.data;fly33.visual_params.data;fly34.visual_params.data;fly35.visual_params.data;...
    fly36.visual_params.data;fly37.visual_params.data;fly38.visual_params.data;fly39.visual_params.data;fly40.visual_params.data;...
    fly41.visual_params.data;fly42.visual_params.data;fly43.visual_params.data;fly44.visual_params.data;fly45.visual_params.data;...
    fly46.visual_params.data;fly47.visual_params.data;fly48.visual_params.data;fly49.visual_params.data;fly50.visual_params.data;...
    fly51.visual_params.data;fly52.visual_params.data;fly53.visual_params.data;fly54.visual_params.data;fly55.visual_params.data;...
    fly56.visual_params.data;fly57.visual_params.data;fly58.visual_params.data;fly59.visual_params.data;fly60.visual_params.data;...
    fly61.visual_params.data;fly62.visual_params.data;fly63.visual_params.data;fly64.visual_params.data;fly65.visual_params.data;...
    fly66.visual_params.data;fly67.visual_params.data;fly68.visual_params.data;fly69.visual_params.data;fly70.visual_params.data;...
    fly71.visual_params.data;fly72.visual_params.data;fly73.visual_params.data;fly74.visual_params.data;fly75.visual_params.data;...
    fly76.visual_params.data;fly77.visual_params.data;fly78.visual_params.data;fly79.visual_params.data;fly80.visual_params.data];

visual.data = alldata;

% Define Fixation Frames
fixation_threshold = 179; % no cut off for fixation
time_thresh = 4; % seconds - did 4 seconds since 170 is roughly 4x 40 fps 
fixation_frames = find(~isnan(visual.data(:,2))); % only calculating on frames where there are no nans
binary_tracking = zeros(length(visual.data),1);
binary_tracking(fixation_frames)=1;
[bouts,lens] = detect_binarybouts_v2(binary_tracking);


% Join "bouts" with <5 frames in between:
flags = zeros(length(bouts),1);
for i = 1:length(bouts)-1
    if i < length(bouts)
        if bouts(i+1, 1)-bouts(i,2) < FPS/5 %any bouts separate by < 200 ms
           binary_tracking(bouts(i,2):bouts(i+1,1)) = ones(length(bouts(i,2):bouts(i+1,1)),1);
        end
    end 
end


% If bouts are less than defined time threshold, exclude:
[bouts,lens] = detect_binarybouts_v2(binary_tracking); % re-calculate bout lengths
for i = 1:length(lens)
    if lens(i)< FPS*time_thresh
        binary_tracking(bouts(i,1):bouts(i,2)) = zeros(length([bouts(i,1):bouts(i,2)]),1);
    end
end

fixation_frames = find(binary_tracking==1); % re-define with processed binary_tracking vector
nbins = round(sqrt(length(fixation_frames)),0);

% Size Histograms -- all assays:

figure(201); clf
nbins = 500;
histogram(visual.data(fixation_frames,4),nbins,'FaceColor', 'k','FaceAlpha',0.2, 'EdgeAlpha',0)
box off; ylabel('Frames','FontSize', 10); xlabel(['Angular Diameter in Y (' char(176) ')'],'FontSize',10); xlim([0 90])
title(['Mean = ' num2str(nanmean(visual.data(fixation_frames,4))) char(176) ' and SD = ' num2str(nanstd(visual.data(fixation_frames,4))) char(176)])
set(gca, 'XScale', 'log'); xlim([0 90]);xticks([2 4.5 9 15 30 45 90]); set(gca, 'XScale', 'log');
for i = 1
    ylim([0 5000])
    if i == 2
        xlim([9 90]) % will line up with the auto setting for the density plot
    end
    if i == 1 
        xlim([9 90])
    end
    set(gca,'LineWidth',2)
    set(gca,'FontSize',14)
end


figure(202); clf
histogram(visual.data(fixation_frames,3),nbins,'FaceColor', 'k','FaceAlpha',0.2, 'EdgeAlpha',0)
box off; ylabel('Frames','FontSize', 10); xlabel(['Angular Diameter in X (' char(176) ')'],'FontSize',10);xlim([0 90])
title(['Mean = ' num2str(nanmean(visual.data(fixation_frames,3))) char(176) ' and SD = ' num2str(nanstd(visual.data(fixation_frames,3))) char(176)])
set(gca, 'XScale', 'log'); xlim([0 90]); xticks([2 4.5 9 15 30 45 90]); set(gca, 'XScale', 'log');
for i = 1
    ylim([0 5000])
    if i == 2
        xlim([9 90]) %will line up with the auto setting for the density plot
    end
    if i == 1 
        xlim([9 90])
    end
    set(gca,'LineWidth',2)
    set(gca,'FontSize',14)
end

% Angular Positions, Velocities, and Looming Velocities in X (Width) and Y (Height):

figure(203); clf
nbins = 500;
subplot(2,1,1)

histogram(visual.data(fixation_frames,2),nbins,'FaceColor', 'k','FaceAlpha',0.2, 'EdgeAlpha',0)
box off; ylabel('Frames','FontSize', 10); xlabel(['Angular Position of Female (' char(176) ')'],'FontSize',14); xlim([-180 180]); ylim([0 1000])
title(['Mean = ' num2str(nanmean(visual.data(fixation_frames,2))) char(176) ' and SD = ' num2str(nanstd(visual.data(fixation_frames,2))) char(176)])

subplot(2,1,2)
histogram(visual.data(fixation_frames,6),nbins,'FaceColor', 'k','FaceAlpha',0.2, 'EdgeAlpha',0)
box off; ylabel('Frames','FontSize', 14); xlabel(['Angular Velocity of Female (' char(176) '/s)'],'FontSize',14); xlim([-300 300])
title(['Mean of Abs Val = ' num2str(nanmean(abs(visual.data(fixation_frames,6)))) char(176) ' and SD = ' num2str(nanstd(abs(visual.data(fixation_frames,6)))) char(176)])
for i =1:2
    subplot(2,1,i)
    box off
    set(gca,'LineWidth',1)
    set(gca,'FontSize',14)
end


figure(204); clf
subplot(2,1,1)
histogram(visual.data(fixation_frames,7),nbins,'FaceColor', 'k','FaceAlpha',0.2, 'EdgeAlpha',0)
box off; ylabel('Frames','FontSize', 14); xlabel(['Angular Expansion of Female Width (' char(176) '/s)'],'FontSize',14); xlim([-300 300])

subplot(2,1,2)
histogram(visual.data(fixation_frames,8),nbins,'FaceColor', 'k','FaceAlpha',0.2, 'EdgeAlpha',0)
box off; ylabel('Frames','FontSize', 14); xlabel(['Angular Expansion of Female Height (' char(176) '/s)'],'FontSize',14); xlim([-300 300])

for i =1:2
    subplot(2,1,i)
    box off
    set(gca,'LineWidth',1)
    set(gca,'FontSize',14)
end



% Size Heatmaps -- all assays:

densityplot(visual.data(fixation_frames,3), visual.data(fixation_frames,4), 'nbins', [30 30])
%colormap(flipud(gray))
set(gca,'Color','k')
colormap(inferno)
xlabel(['Width (' char(176) ')'], 'FontSize', 14);
ylabel(['Height (' char(176) ')'], 'FontSize', 14);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xticks([2 4.5 9 15 30 45 90])
yticks([2 4.5 9 15 30 45 90])
set(gca,'TickDir','out');
axis equal

xlim([0 90])
ylim([0 90])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca,'LineWidth',2)
set(gca,'FontSize',14)
%axis equal



% Reference Figure: 
% plot curves that illustrate all possible female angular widths as a function of
% inter-fly distance (if approximate a female as an ellipse with major axis
% = 3 mm and minor axis = 1 mm). In this sense, the lower bound is also an
% an approximation of all possible angular heights, if approximate the
% female's height to be 1 mm

dists = 0:0.1:40;
minWidths = 2.*atand(0.5./dists);
maxWidths = 2.*atand(1.5./dists);

figure
patch([dists fliplr(dists)], [minWidths fliplr(maxWidths)], 'k','FaceAlpha', 0.2); hold on
plot(dists,minWidths,'c','LineWidth',2); hold on
plot(dists,maxWidths,'m','LineWidth',2); hold on
ylabel(['Angular Size (' char(176) ')'], 'FontSize', 14);
set(gca,'LineWidth',2)
set(gca,'FontSize',14)
xlabel('Inter-fly Distance (mm)')
legend('All Possible Widths','Height or Minimum Width','Maximum Width')
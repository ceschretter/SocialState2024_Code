
%% Descriptive Analysis of Visual Experience of Male during Courtship
% Load workspace with calibration,tracks, and features from FlyTracker:

%load('melanogaster_flytrackeroutput.mat');
load('melanogaster_courtshipdyads_flytrackeroutput.mat');

%% Generate Visual Feature in pre-copulation frames for each assay at 40 Hz Sampling Rate:

% collect visual features from courtship start to copulation (or end of
% assay):

nose = false; %set to true if want to calculate based on anterior-most point of fly; set to false if want to calculate
                % visual features W.R.T. to subject's centroid:
P1_visfeats = GetVisualFeatures(P1_tracks,P1_calib,1,2,40,[P1_calib.start P1_calib.cop_ind],nose,false); 
P2_visfeats = GetVisualFeatures(P2_tracks,P2_calib,1,2,40,[P2_calib.start P2_calib.cop_ind],nose,false);
P3_visfeats = GetVisualFeatures(P3_tracks,P3_calib,1,2,40,[P3_calib.start P3_calib.cop_ind],nose,false);
P4_visfeats = GetVisualFeatures(P4_tracks,P4_calib,1,2,40,[P4_calib.start length(P4_tracks.data)],nose,false); 
P5_visfeats = GetVisualFeatures(P5_tracks,P5_calib,1,2,40,[P5_calib.start P5_calib.cop_ind],nose,false);
P6_visfeats = GetVisualFeatures(P6_tracks,P6_calib,1,2,40,[P6_calib.start P6_calib.cop_ind],nose,false);
P7_visfeats = GetVisualFeatures(P7_tracks,P7_calib,1,2,40,[P7_calib.start length(P7_tracks.data)],nose,false);
P8_visfeats = GetVisualFeatures(P8_tracks,P8_calib,1,2,40,[P8_calib.start P8_calib.cop_ind],nose,false);
P9_visfeats = GetVisualFeatures(P9_tracks,P9_calib,1,2,40,[P9_calib.start P9_calib.cop_ind],nose,false);
P10_visfeats = GetVisualFeatures(P10_tracks,P10_calib,1,2,40,[P10_calib.start P10_calib.cop_ind],nose,false);
P11_visfeats = GetVisualFeatures(P11_tracks,P11_calib,1,2,40,[P11_calib.start P11_calib.cop_ind],nose,false);
P12_visfeats = GetVisualFeatures(P12_tracks,P12_calib,1,2,40,[P12_calib.start P12_calib.cop_ind],nose,false);
P13_visfeats = GetVisualFeatures(P13_tracks,P13_calib,1,2,40,[P13_calib.start P13_calib.cop_ind],nose,true);

%% Plot Visual Features across all flies:
% Plot female size distributions from 13 assays:

all_data = [P1_visfeats.data; P2_visfeats.data; P3_visfeats.data; P4_visfeats.data; P5_visfeats.data;...
    P6_visfeats.data; P7_visfeats.data; P8_visfeats.data; P9_visfeats.data; P10_visfeats.data; P11_visfeats.data;...
    P12_visfeats.data; P13_visfeats.data];

visual.data = all_data; 

%% If want to define fixation frames within certain range of visual field or for minimum time:
fixation_threshold = 179; % cut-off for fixation is +/- 15 deg relative to male's visual field center (0 deg)
time_thresh = 0; % seconds
fixation_frames = find(abs(visual.data(:,2))< fixation_threshold); % if want to look at experience when NOT fixated, just change '<' to '>'
binary_tracking = zeros(length(visual.data),1);
binary_tracking(fixation_frames)=1;
[bouts,lens] = detect_binarybouts_v2(binary_tracking);
FPS = 40;

% Join "bouts" with <5 frames in between:
flags = zeros(length(bouts),1);
for i = 1:length(bouts)-1
    if i < length(bouts)
        if bouts(i+1, 1)-bouts(i,2) < FPS/5 %any bouts separate by < 200 ms
           %bouts(i,2) = bouts(i+1,2);
           %bouts(i+1,:) = [];
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

%% Size Histograms -- all assays:

figure(201); clf
subplot(2,1,1)
%nbins = 100;
histogram(visual.data(fixation_frames,4),nbins,'FaceColor', 'm','FaceAlpha',1, 'EdgeAlpha',0)
box off; ylabel('Frames','FontSize', 10); xlabel(['Angular Diameter in Y (' char(176) ')'],'FontSize',10); xlim([0 90])
title(['Mean = ' num2str(mean(visual.data(fixation_frames,4))) char(176) ' and SD = ' num2str(std(visual.data(fixation_frames,4))) char(176)])
set(gca, 'XScale', 'log'); xlim([0 90]);xticks([2 4.5 9 15 30 45 90]); set(gca, 'XScale', 'log');
subplot(2,1,2)
histogram(visual.data(fixation_frames,3),nbins,'FaceColor', 'm','FaceAlpha',1, 'EdgeAlpha',0)
box off; ylabel('Frames','FontSize', 10); xlabel(['Angular Diameter in X (' char(176) ')'],'FontSize',10);xlim([0 90])
title(['Mean = ' num2str(mean(visual.data(fixation_frames,3))) char(176) ' and SD = ' num2str(std(visual.data(fixation_frames,3))) char(176)])
set(gca, 'XScale', 'log'); xlim([0 90]); xticks([2 4.5 9 15 30 45 90]); set(gca, 'XScale', 'log');
for i = 1:2
    subplot(2,1,i)
    ylim([0 10000])
    if i == 2
        xlim([2 90]) % the density plot is stuck at auto settings for xlim and ylim, so this will line up with that auto setting
    end
    if i == 1 
        xlim([2 90])
    end
    set(gca,'LineWidth',2)
    set(gca,'FontSize',14)
end

%% Angular Positions, Velocities, and Looming Velocities in X (Width) and Y (Height):

figure(203); clf
subplot(4,1,1)

histogram(visual.data(fixation_frames,2),nbins,'FaceColor', 'm','FaceAlpha',1, 'EdgeAlpha',0)
box off; ylabel('Frames','FontSize', 10); xlabel(['Angular Position of Female (' char(176) ')'],'FontSize',14); xlim([-180 180])
title(['Mean = ' num2str(mean(visual.data(fixation_frames,2))) char(176) ' and SD = ' num2str(std(visual.data(fixation_frames,2))) char(176)])

subplot(4,1,2)
histogram(visual.data(fixation_frames,6),nbins,'FaceColor', 'm','FaceAlpha',1, 'EdgeAlpha',0)
box off; ylabel('Frames','FontSize', 14); xlabel(['Angular Velocity of Female (' char(176) '/s)'],'FontSize',14); xlim([-300 300])
title(['Mean of Abs Val = ' num2str(mean(abs(visual.data(fixation_frames,6)))) char(176) ' and SD = ' num2str(std(abs(visual.data(fixation_frames,6)))) char(176)])

subplot(4,1,3)
histogram(visual.data(fixation_frames,7),nbins,'FaceColor', 'm','FaceAlpha',1, 'EdgeAlpha',0)
box off; ylabel('Frames','FontSize', 14); xlabel(['Angular Expansion of Female Width (' char(176) '/s)'],'FontSize',14); xlim([-300 300])

subplot(4,1,4)
histogram(visual.data(fixation_frames,8),nbins,'FaceColor', 'm','FaceAlpha',1, 'EdgeAlpha',0)
box off; ylabel('Frames','FontSize', 14); xlabel(['Angular Expansion of Female Height (' char(176) '/s)'],'FontSize',14); xlim([-300 300])

for i =1:4
    subplot(4,1,i)
    box off
    set(gca,'LineWidth',1)
    set(gca,'FontSize',14)
end




%% Size Heatmaps -- all assays:

densityplot(visual.data(fixation_frames,3), visual.data(fixation_frames,4), 'nbins', [120  120])
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
%clim([0 1500])
xlim([2 90])
ylim([2 90])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca,'LineWidth',2)
set(gca,'FontSize',14)
%axis equal



%% Reference Figure: 
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









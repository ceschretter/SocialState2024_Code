function [xe,ye] = get_ellipse(trk, ctr_coords, theta, fly_ind)
% Input arguments
%    prefix = file prefix for -track.mat (usually the same prefix for both files)

%    ctr_f = x,y coordinates of female; typically nx2 matrix, for n frames

%    theta = orientation of female

%    fly_ind: integer corresponding to FlyTracker ID for fly of interest


% Output arguments

%    xe = x coordinates of points around perimeter of ellipse fit to female
%    ye = y coordinates of points around perimeter of ellipse fit to female

%% load track.mat and extract parameters:

trk_struc = trk;

x = ctr_coords(:,1); % fly x position
y = ctr_coords(:,2); % fly y position
major_axis_len = trk_struc.data(fly_ind,:,4); % fly major axis length
minor_axis_len = trk_struc.data(fly_ind,:,5); % fly minor axis length
a = major_axis_len./2;
b = minor_axis_len./2;

%% fit ellipse, extract perimeter coordinates! 
%  Adapted from https://www.mathworks.com/matlabcentral/answers/129273-to-plot-an-ellipse

clear phi X1 Y1
phi = 0:0.1:2*pi;
costheta = cos(theta);
sintheta = sin(theta);

xe = zeros(length(theta), length(phi)); ye = zeros(length(theta), length(phi));
for i = 1:length(theta)
    X1 = a(1,i)*cos(phi);
    Y1 = b(1,i)*sin(phi);
    rot_mat = [costheta(i) -sintheta(i); sintheta(i) costheta(i)];
    for j = 1:length(phi)
        rot_ellipse(j,:) = rot_mat*([X1(1,j) Y1(1,j)]');
    end
    xe(i,:) = rot_ellipse(:,1)' + x(i);
    ye(i,:) = rot_ellipse(:,2)' + y(i);
end

display('Target fly elliptical fits complete for every frame!');
% 
%% Plots -- can use this to plot if desired
% if plotit > 0 % plot ellipse if desired 
%         figure
%         for i = 1:end_ind
%             plot(xe(i,:),ye(i,:),'-','Color',[0 0 1],'LineWidth',2);
%             xlim([-1500 1500]); ylim([-1500 1500])
%             line([-1500 1500], [0 0], 'Color', 'k'); line([0 0], [-1500 1500],  'Color', 'k')
%         end
% end

end


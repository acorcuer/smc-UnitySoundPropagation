clear all;
impulse = [1 zeros(1, 99)];
fs = 44100;
c = 340;
theta = 45;
% Position of the sound: x, y (only in 2D for the moment) in meters
x0 = 2;
y0 = 2;
% Listener
dx_l = 3;
dy_l = 3;
listener = [dx_l dy_l];
scatter(x0, y0, 45, 'filled');hold on
scatter(dx_l, dy_l, 45, 'filled')
radius_list = 2;
% Direct sound
    dx = abs(x0 - dx_l);
    dy = abs(y0 - dy_l);
    distance = sqrt(dx*dx + dy*dy);
    delay = distance/c;
    m_direct = 1/distance;
    d_direct = floor(delay * fs);
% Number of reflections: n_reflections
n_reflections = 200;
attenuation = 0.5;
h = zeros(1,n_reflections);
att_ray = 1;
delays = zeros(1, n_reflections+1);
delays(1) = d_direct;
k=2;
h(delays(1)) = m_direct;
mags = zeros(1, n_reflections);
% General equation of a line Ax+By+C=0
j=1;
walls = [0 1 -4; 1 0 -1; 0 1 -1; 1 0 -6];
l2 = [-(sind(theta))  (cosd(theta))  ((-cosd(theta))*y0 + (sind(theta))*x0)];
m = -l2(1)/l2(2); b = -l2(3)/l2(2); x = 1:10;
            plot(x, m*x+b); hold on
for i=1:(n_reflections+1)
    while j < 5
        l1 = walls(j, :);
        % Intersection point
        intersection_x = det([l1(2:3); l2(2:3)])/det([l1(1:2); l2(1:2)]);
        intersection_y = det([l1(3) l1(1); l2(3) l2(1)])/det([l1(1:2); l2(1:2)]);
        % Jordan curve theorem: TODO
        if intersection_x <= -walls(4,3) && intersection_y <= -walls(1,3) && ...
                intersection_y >= -walls(3,3) && intersection_x >= -walls(2,3)...
                && intersection_x ~= x0 && intersection_y ~= y0
            a = [intersection_x intersection_y] - [x0 y0];
            b = listener - [x0 y0];
            dist_pointLine = abs( det([a-b;c-b]) )/norm(c-b) ;
            % Length of the ray
            dx = abs(x0-intersection_x);
            dy = abs(y0-intersection_y);
            dist = sqrt(dx*dx + dy*dy);
            % New source
            scatter(intersection_x,intersection_y); hold on
            x0 = intersection_x;
            y0 = intersection_y;
            % Reflection
           % r = [-l2(2) l2(1)]+2*cosd(theta)*[l1(1) l1(2)];
           N = [l1(1) l1(2)];
           N = N / norm(N);
           r = [-l2(2) l2(1)] - 2*dot([-l2(2) l2(1)],N)*N;
            l2=[r(2) -r(1) (r(1)*y0 - r(2)*x0)];
            m = -l2(1)/l2(2); b = -l2(3)/l2(2); x = 1:10;
            plot(x, m*x+b);hold on
            j = j + 5;
        end
        j = j+1;
    end
    % Let's check if the ray has crossed the receiver
    if dist_pointLine < radius_list
        delays(k) = ceil((dist/c) * fs);
        att_ray = att_ray * attenuation;
        mags(k) = att_ray;
        h(delays(k)) = att_ray;
        k=k+1;
    end
    j=1;
end

% figure
% stem(h)

%     N = N / norm(N);
%     m = -l2(1)/l2(2); b = -l2(3)/l2(2); x = 1:10;
%     plot(x, m*x+b);hold on
% 
% scatter(r(1),r(2))
% scatter(intersection_x,intersection_y)
plot(1:6,4*ones(1,6), 'LineWidth', 2);hold on
plot(1:6,1*ones(1,6), 'LineWidth', 2);hold on
plot(1*ones(1,4), 1:4, 'LineWidth', 2);hold on
plot(6*ones(1,4), 1:4, 'LineWidth', 2)
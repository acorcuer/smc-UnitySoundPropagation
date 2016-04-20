impulse = [1 zeros(1, 99)];
fs = 44100;
c = 340;
theta = 30;
% Position of the sound: x, y (only in 2D for the moment) in meters
x0 = 2;
y0 = 2;
% Listener
dx_l = 3;
dy_l = 2;
listener = [dx_l dy_l];
radius_list = 4;
% Direct sound
    dx = abs(x0 - dx_l);
    dy = abs(y0 - dy_l);
    distance = sqrt(dx*dx + dy*dy);
    delay = distance/c;
    m_direct = 1/distance;
    d_direct = floor(delay * fs);
% Number of reflections: n_reflections
n_reflections = 5000;
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
for i=1:(n_reflections+1)
    l2 = [(sind(theta) - y0)  -(cosd(theta) - x0)  ((x0-cosd(theta))*y0 + (sind(theta)-y0)*x0)];
    while j < 5
        l1 = walls(j, :);
        % Intersection point
        intersection_x = det([l1(2:3); l2(2:3)])/det([l1(1:2); l2(1:2)]);
        intersection_y = det([l1(3) l1(1); l2(3) l2(1)])/det([l1(1:2); l2(1:2)]);
        % Jordan curve theorem: TODO
        if intersection_x <= -walls(4,3) && intersection_y <= -walls(1,3) && ...
                intersection_y >= -walls(3,3) && intersection_x >= -walls(2,3)
            a = [intersection_x intersection_y] - [x0 y0];
            b = listener - [x0 y0];
            dist_pointLine = abs( det([a-b;c-b]) )/norm(c-b) ;
            % Length of the ray
            dx = abs(x0-intersection_x);
            dy = abs(y0-intersection_y);
            dist = sqrt(dx*dx + dy*dy);
            % New source
            x0 = intersection_x;
            y0 = intersection_y;
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

stem(h)
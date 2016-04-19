% [audio, fs] = audioread('Piano.wav');
impulse = [1 zeros(1, 99)];
fs = 44100;
distance = 20;
c = 340;

delay = distance/c;
gain = 1/distance;
source = impulse*gain;
d_direct = floor(delay * fs);
theta = 30;
% Position of the sound: x, y (only in 2D for the moment) in meters
x0 = 2;
y0 = 2;
% Listener
dx_l = 3;
dy_l = 2;
listener = [dx_l dy_l];
radius_list = 2;
% Number of reflections: n_reflections
n_reflections = 10;
attenuation = 0.5;
h = zeros(1,10000);
att_ray = 1;
delays = zeros(1, n_reflections+1);
delays(1) = d_direct;
mags = zeros(1, n_reflections);
% General equation of a line Ax+By+C=0
walls = [0 1 4; 1 0 1; 0 1 1; 1 0 6];
for i=1:i<n_reflections + 1
%     dx = abs(dx_s - dx_s1);
%     dy = abs(dy_s - dy_s1);
    k=1;
    j = 1;
    l2 = [(sind(theta) - y0)  (cosd(theta) - x0)  ((x0-cosd(theta))*y0 + (sind(theta)-y0)*x0)];
    while j < 5
        l1 = walls(j, :);
        % Intersection point
        intersection_x = det([l1(2:3); l2(2:3)])/det([l1(1:2); l2(1:2)]);
        intersection_y = det([l1(3) l1(1); l2(3) l2(1)])/det([l1(1:2); l2(1:2)]);
        % Jordan curve theorem: TODO
        if intersection_x < walls(4,3) && intersection_y < walls(1,3) && ...
                intersection_y > walls(3,3) && intersection_x > walls(2,3)
            a = [intersection_x intersection_y] - [x0 y0];
            b = listener - [x0 y0];
            dist_pointLine = norm(cross(a,b)) / norm(a);
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
        delays(k) = floor(dist/c) * fs;
        att_ray = att_ray * attenuation;
        mags(k) = att_ray;
        h(delays(k)) = att_ray;
        k=k+1;
    end

end

   
    
    
    
    
    

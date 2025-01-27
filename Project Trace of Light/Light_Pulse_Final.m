% Initialization
clear
clc
close all

% Definition of the range of our render
N = 1000;
x = -10:20/N:10-0.04;
y = -10:20/N:10-0.04;
f = 150; % focale
r = 1.5; % pupil
distance = -20:1:100; % Propagation Distance

% Information of a light pulse
t_propagation = 0.3e-11;
t_duration = 0.1e-11;
t_gap = 0.2e-11;

% Three different colors
lambda_r = 0.65e-3;
lambda_g = 0.55e-3;
lambda_b = 0.45e-3;

% Calculate the propagation for three different colors
propagation_r = Propagation_in_free_space(lambda_r,distance,x,y,f,r,t_propagation,t_duration,t_gap);
propagation_g = Propagation_in_free_space(lambda_g,distance,x,y,f,r,t_propagation,t_duration,t_gap);
propagation_b = Propagation_in_free_space(lambda_b,distance,x,y,f,r,t_propagation,t_duration,t_gap);

% 3D graph
Draw_3D_graph(x,y,distance,propagation_r,'r');
Draw_3D_graph(x,y,distance,propagation_g,'g');
Draw_3D_graph(x,y,distance,propagation_b,'b');
Draw_3D_graph(x,y,distance,(propagation_r+propagation_g+propagation_b)/3,'c');

% 2D slice graph along the light axis
Draw_slice_plot(propagation_r,'r');
Draw_slice_plot(propagation_g,'g');
Draw_slice_plot(propagation_b,'b');
Draw_slice_plot((propagation_r+propagation_g+propagation_b)/3,'c')

% 1D intensity graph along the light axis
figure;
hold on
Draw_intensity_plot(propagation_r,distance,lambda_r,'r');
Draw_intensity_plot(propagation_g,distance,lambda_g,'g');
Draw_intensity_plot(propagation_b,distance,lambda_b,'b');
hold off

% In this function, we calculate the light propagation in the free space
% First we do the FFt to get the Airy disk in the focal plane
% Second we do anthor FFT to tranfer the graph in position domain to frequency domain
% Then we multiply the transfer function to calculate the light distribution in different distance (in frequency domain)
% Lastly we do inverse FFT to get the intensity distribution in position  domain
function propagation = Propagation_in_free_space(lambda,distance,x,y,f,r,t_propagation,t_duration,t_gap)
    propagation = []; % save the intensity in different domain
    c = 3 * 10e11; % light speed in mm
    light_range = -20 + t_propagation*c : 1 : -20 + t_duration*c + t_propagation*c; % light pluse range
    distance_gap = t_gap*c; % distance bewteen two pluse
    max_n = 120 / distance_gap;  % amount of light pluse in our range
    i = 0;
    k = 2 * pi / lambda; % light vector

    [X,Y] = meshgrid(x,y);
    P1 = (X.^2 + Y.^2) < r^2;

    % Aberration Spherique
    r_grid = sqrt(X.^2 + Y.^2);
    alpha = 2;
    phi_spherical = alpha * r_grid.^4;

    % Coma calculation
    r_grid(r_grid == 0) = eps;
    cos_theta = X ./ r_grid;
    gamma = 2;
    phi_coma = gamma * r_grid.^3 .* cos_theta;

    % Field Curvature
    beta = 0.5;
    phi_field_curvature = beta * r_grid.^2;

    P1_with_lens = P1 *exp(1j * (phi_spherical + phi_coma + phi_field_curvature) ) ;
    Airy = fftshift(fft2(P1_with_lens));
    P2 = fft2(Airy);

    for d = distance
        vx = x / (lambda * f);
        vy = y / (lambda * f);
        [Vx, Vy] = meshgrid(vx, vy);
    
        is_in_range = false;
        for n = 0:max_n
        % Draw all the pulse
            shifted_min = min(light_range) + n * distance_gap;
            shifted_max = max(light_range) + n * distance_gap;
    
            if d >= shifted_min && d <= shifted_max
                is_in_range = true;
                break;
            end
        end
    
        if is_in_range
            if d == 0
                P3 = Airy;
            else
                Hd = exp(-1j * k * d) .* exp(1j * pi * lambda * d * (Vx.^2 + Vy.^2));
                P3 = ifft2(P2 .* Hd);
            end
        else
            P3 = zeros(size(Airy));
        end
    
        i = i + 1;
        propagation(:, :, i) = P3; % propagation result
    end
end

% This function draw the 2D slice of the result
function Draw_slice_plot(propagation,color)
    y = -10:20/1000:10-0.04;
    X_slice = squeeze(propagation(500, :, :));
    Y_slice = squeeze(propagation(:, 500, :));
    x_slice_slice = X_slice(501, :);
    y_slice_slice = Y_slice(501, :);
    [m,n] = size(X_slice);

    figure();
    VisuIdBPh([],y,X_slice,-40);
    xlabel('Distance Along Z axis');
    ylabel('Pupille');
    title(sprintf('Propagation for the Light, Color: %s', color));

end

% This function draw the Intensity vs Optical Axis graph
function Draw_intensity_plot(propagation,distance,lambda,color)
    X_slice = squeeze(propagation(500, :, :));
    Y_slice = squeeze(propagation(:, 500, :));
    x_slice_slice = X_slice(501, :);
    y_slice_slice = Y_slice(501, :);
    [X_max, linear_idx] = max(abs(X_slice(:)));
    [row_idx, col_idx] = ind2sub(size(X_slice), linear_idx);
    max_distance = distance(col_idx);
    disp(['Maxmimum Column：', num2str(row_idx)]);
    disp(['Maximum Line：', num2str(col_idx)]);
    disp(['Maximum Value：', num2str(X_max)]);

    plot(distance, abs(x_slice_slice), color, 'DisplayName', ['\lambda = ', num2str(lambda), ' nm']);
    xlabel('Distance along Z axis')
    ylabel('Intensity in the center of Pupille')
    title(sprintf('Intensity with propagation distance, Color: %s', color))
end

% This function draw the slice in the 3D graph
function Draw_3D_graph(x,y,distance,propagation,color)
    [X,Y,Z] = meshgrid(x, y, distance);
    figure('position',[100 50 560 780],'color','w');
    h = slice(X,Y,Z,abs(propagation),0,0,[0 40 80]);
    set(h,'facealpha',0.5); 
    shading interp;
    view(-55,5);
    colormap(jet);
    title(sprintf('Propagation for the Light, Color: %s', color), 'FontSize', 14, 'Color', 'k');
end



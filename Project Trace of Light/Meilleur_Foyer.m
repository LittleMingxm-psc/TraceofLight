clear
clc
close all

N = 1000;
x = -10:20/N:10-0.04;
y = -10:20/N:10-0.04;
f = 150;

r = 1.5;
distance = -40:1:100;

lambda_r = 0.65e-3;
lambda_g = 0.55e-3;
lambda_b = 0.45e-3;

propagation_r = Propagation_in_free_space(lambda_r,distance,x,y,f,r);
propagation_g = Propagation_in_free_space(lambda_g,distance,x,y,f,r);
propagation_b = Propagation_in_free_space(lambda_b,distance,x,y,f,r);

propagation_r_meilleur = meilleur_foyer(lambda_r,distance,x,y,f,r);
propagation_g_meilleur = meilleur_foyer(lambda_g,distance,x,y,f,r);
propagation_b_meilleur = meilleur_foyer(lambda_b,distance,x,y,f,r);

Draw_3D_graph(x,y,distance,propagation_r,'r');
Draw_3D_graph(x,y,distance,propagation_g,'g');
Draw_3D_graph(x,y,distance,propagation_b,'b');
Draw_3D_graph(x,y,distance,(propagation_r+propagation_g+propagation_b)/3,'c');

Draw_slice_plot(propagation_r,'r');
Draw_slice_plot(propagation_g,'g');
Draw_slice_plot(propagation_b,'b');
Draw_slice_plot((propagation_r+propagation_g+propagation_b)/3,'c')

figure;
hold on
Draw_intensity_plot(propagation_r,distance,lambda_r,'r');
Draw_intensity_plot(propagation_g,distance,lambda_g,'g');
Draw_intensity_plot(propagation_b,distance,lambda_b,'b');
hold off

% 1D intensity graph along the light axis after the meilleur_foyer
figure;
hold on
Draw_intensity_plot(propagation_r_meilleur,distance,lambda_r,'r');
Draw_intensity_plot(propagation_g_meilleur,distance,lambda_g,'g');
Draw_intensity_plot(propagation_b_meilleur,distance,lambda_b,'b');
hold off

function propagation = Propagation_in_free_space(lambda,distance,x,y,f,r)
    propagation = [];
    i = 0;
    k = 2 * pi / lambda;

   [X,Y] = meshgrid(x,y);
    P1 = (X.^2 + Y.^2) < r^2;

    % Aberration Spherique
    r_grid = sqrt(X.^2 + Y.^2);
    alpha = 5;
    phi_spherical = alpha * r_grid.^4;

    % Coma calculation
    r_grid(r_grid == 0) = eps;
    cos_theta = X ./ r_grid;
    gamma = 0.1;
    phi_coma = gamma * r_grid.^3 .* cos_theta;

    % Field Curvature
    beta = 0.2;
    phi_field_curvature = beta * r_grid.^2;

    P1_with_lens = P1 .* exp(1j *(phi_spherical + phi_coma + phi_field_curvature));
    Airy = fftshift(fft2(P1_with_lens)); 

    P2 = fft2(Airy);
    imagesc(abs(P2))

    for d = distance
        vx = x / (lambda * f);
        vy = y / (lambda * f);
        [Vx, Vy] = meshgrid(vx, vy);
        if d == 0
            P3 = Airy;
        else
            Hd = exp(-1j * k * d) .* exp(1j * pi * lambda * d * (Vx.^2 + Vy.^2));
            P3 = ifft2(P2 .* Hd);
        end
        i = i + 1;
        propagation(:, :, i) = P3;
    end
end

% Calulate for the meilleur foyer
function propagation = meilleur_foyer(lambda,distance,x,y,f,r)
    propagation = [];
    i = 0;
    k = 2 * pi / lambda;

   [X,Y] = meshgrid(x,y);
    P1 = (X.^2 + Y.^2) < r^2;

    % The formula for the Aberrationa are changed to meilleur foyer one
    % Aberration Spherique
    r_grid = sqrt(X.^2 + Y.^2);
    alpha = 5;
    phi_spherical = alpha *( r_grid.^4 - r_grid.^2);

    % Coma calculation
    r_grid(r_grid == 0) = eps;
    cos_theta = X ./ r_grid;
    gamma = 0.1;
    phi_coma = gamma * (r_grid.^3 - 2/3*r_grid) .* cos_theta;

    % Field Curvature
    beta = 0.2;
    phi_field_curvature = beta * r_grid.^2;

    P1_with_lens = P1 .* exp(1j *(phi_spherical + phi_coma + phi_field_curvature));
    Airy = fftshift(fft2(P1_with_lens));

    P2 = fft2(Airy);
    imagesc(abs(P2))

    for d = distance
        vx = x / (lambda * f);
        vy = y / (lambda * f);
        [Vx, Vy] = meshgrid(vx, vy);
        if d == 0
            P3 = Airy;
        else
            Hd = exp(-1j * k * d) .* exp(1j * pi * lambda * d * (Vx.^2 + Vy.^2));
            P3 = ifft2(P2 .* Hd);
        end
        i = i + 1;
        propagation(:, :, i) = P3; 
    end
end

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

function Draw_intensity_plot(propagation,distance,lambda,color)
    X_slice = squeeze(propagation(500, :, :));
    Y_slice = squeeze(propagation(:, 500, :));
    x_slice_slice = X_slice(501, :);
    y_slice_slice = Y_slice(501, :);
    [X_max, linear_idx] = max(abs(X_slice(:)));
    [row_idx, col_idx] = ind2sub(size(X_slice), linear_idx);
    max_distance = distance(col_idx);
    disp(['Maximum column：', num2str(row_idx)]);
    disp(['Maximum line：', num2str(col_idx)]);
    disp(['Maximum Intensity：', num2str(X_max)]);

    plot(distance, abs(x_slice_slice), color, 'DisplayName', ['\lambda = ', num2str(lambda), ' nm']);
    xlabel('Distance along Z axis')
    ylabel('Intensity in the center of Pupille')
    title(sprintf('Intensity with propagation distance'))
end

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



% Add the path to the Python script
insert(py.sys.path, int32(0), '.');  % Assuming the Python script is in the current directory

% Load data
A = dlmread('log_files/vertex.txt');
N = A(1, 1);
Lx = A(2, 1);
Ly = A(3, 2);
vertex = A(4:3+N, :);
H1 = dlmread('log_files/connectivity_matrix.txt');

% Initialize arrays to store bond lengths and bond angles
bond_lengths = [];
bond_angles = [];

% Plot the lattice and calculate bond lengths
figure;
hold on;
for i = 1:N
    for j = 1:i
        if H1(i, j) == 1
            dx = abs(vertex(i, 1) - vertex(j, 1));
            dy = abs(vertex(i, 2) - vertex(j, 2));
            bond_length = sqrt(dx^2 + dy^2);
            bond_lengths = [bond_lengths; bond_length];
            
            % Apply periodic boundary conditions for plotting
            if (dx < Lx / 2.0) && (dy < Ly / 2.0)
                plot([vertex(i, 1) vertex(j, 1)], [vertex(i, 2) vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
            elseif (dx < Lx / 2.0) && (dy > Ly / 2.0)
                if vertex(i, 2) < vertex(j, 2)
                    plot([vertex(i, 1) vertex(j, 1)], [vertex(i, 2) + Ly vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
                    plot([vertex(i, 1) vertex(j, 1)], [vertex(i, 2) vertex(j, 2) - Ly], 'Color', 'r', 'LineWidth', 2);
                else
                    plot([vertex(i, 1) vertex(j, 1)], [vertex(i, 2) vertex(j, 2) + Ly], 'Color', 'r', 'LineWidth', 2);
                    plot([vertex(i, 1) vertex(j, 1)], [vertex(i, 2) - Ly vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
                end
            elseif (dx > Lx / 2.0) && (dy < Ly / 2.0)
                if vertex(i, 1) < vertex(j, 1)
                    plot([vertex(i, 1) + Lx vertex(j, 1)], [vertex(i, 2) vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
                    plot([vertex(i, 1) vertex(j, 1) - Lx], [vertex(i, 2) vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
                else
                    plot([vertex(i, 1) - Lx vertex(j, 1)], [vertex(i, 2) vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
                    plot([vertex(i, 1) vertex(j, 1) + Lx], [vertex(i, 2) vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
                end
            else
                if (vertex(i, 1) - Lx / 2.0) * (vertex(i, 2) - Ly / 2.0) > 0.0
                    if vertex(i, 1) > vertex(j, 1)
                        plot([vertex(i, 1) vertex(j, 1) + Lx], [vertex(i, 2) vertex(j, 2) + Ly], 'Color', 'r', 'LineWidth', 2);
                        plot([vertex(i, 1) vertex(j, 1) + Lx], [vertex(i, 2) - Ly vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
                        plot([vertex(i, 1) - Lx vertex(j, 1)], [vertex(i, 2) vertex(j, 2) + Ly], 'Color', 'r', 'LineWidth', 2);
                        plot([vertex(i, 1) - Lx vertex(j, 1)], [vertex(i, 2) - Ly vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
                    else
                        plot([vertex(i, 1) + Lx vertex(j, 1)], [vertex(i, 2) + Ly vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
                        plot([vertex(i, 1) vertex(j, 1) - Lx], [vertex(i, 2) + Ly vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
                        plot([vertex(i, 1) + Lx vertex(j, 1)], [vertex(i, 2) vertex(j, 2) - Ly], 'Color', 'r', 'LineWidth', 2);
                        plot([vertex(i, 1) vertex(j, 1) - Lx], [vertex(i, 2) vertex(j, 2) - Ly], 'Color', 'r', 'LineWidth', 2);
                    end
                else
                    if vertex(i, 1) > vertex(j, 1)
                        plot([vertex(i, 1) vertex(j, 1) + Lx], [vertex(i, 2) vertex(j, 2) - Ly], 'Color', 'r', 'LineWidth', 2);
                        plot([vertex(i, 1) - Lx vertex(j, 1)], [vertex(i, 2) vertex(j, 2) - Ly], 'Color', 'r', 'LineWidth', 2);
                        plot([vertex(i, 1) - Lx vertex(j, 1)], [vertex(i, 2) + Ly vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
                        plot([vertex(i, 1) vertex(j, 1) + Lx], [vertex(i, 2) + Ly vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
                    else
                        plot([vertex(i, 1) + Lx vertex(j, 1)], [vertex(i, 2) - Ly vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
                        plot([vertex(i, 1) vertex(j, 1) - Lx], [vertex(i, 2) vertex(j, 2) + Ly], 'Color', 'r', 'LineWidth', 2);
                        plot([vertex(i, 1) vertex(j, 1) - Lx], [vertex(i, 2) - Ly vertex(j, 2)], 'Color', 'r', 'LineWidth', 2);
                        plot([vertex(i, 1) + Lx vertex(j, 1)], [vertex(i, 2) vertex(j, 2) + Ly], 'Color', 'r', 'LineWidth', 2);
                    end
                end
            end
        end
    end
end


% Plot settings
axis equal;
axis([0 Lx 0 Ly]);
set(gca, 'xtick', [], 'ytick', []);
set(gca, 'FontSize', 15);
plot([0.0 0.0], [0.0 Ly], 'Color', 'k', 'LineWidth', 2);
plot([Lx Lx], [0.0 Ly], 'Color', 'k', 'LineWidth', 2);
plot([0.0 Lx], [0.0 0.0], 'Color', 'k', 'LineWidth', 2);
plot([0.0 Lx], [Ly Ly], 'Color', 'k', 'LineWidth', 2);
hold off;

saveas(gcf, 'log_files_plot\structure_plot.png');

% Calculate bond angles
for i = 1:N
    neighbors = find(H1(i, :));
    for j = 1:length(neighbors)
        for k = j+1:length(neighbors)
            v1 = vertex(neighbors(j), :) - vertex(i, :);
            v2 = vertex(neighbors(k), :) - vertex(i, :);
            angle = atan2d(norm(cross([v1 0], [v2 0])), dot(v1, v2));
            bond_angles = [bond_angles; angle];
        end
    end
end

% Plot histograms of bond lengths and bond angles
figure;
subplot(1, 2, 1);
histogram(bond_lengths, 'Normalization', 'count', 'BinWidth', 0.05);
title('Bond Lengths');
xlabel('Bond length (Å)');
ylabel('Counts');
xlim([0.9 2.1]);
ylim([0 300]);

subplot(1, 2, 2);
histogram(bond_angles, 'Normalization', 'count', 'BinWidth', 1);
title('Bond Angles');
xlabel('Bond angle (degrees)');
ylabel('Counts');
xlim([80 160]);
ylim([0 250]);

saveas(gcf, 'log_files_plot\bond_statistics_plot.png');

% Construct adjacency matrix and find rings using Python code
max_ring_size = 9; % Set maximum ring size
cutoff = 2; % Adjust the cutoff distance based on your system

% Convert vertex data to Python format
coords_py = py.numpy.array(vertex);
num_atoms_py = int32(N); % Convert to Python int type

% Call the Python function to find rings
py_script = 'ring_statistics';
py.importlib.import_module(py_script); % Import the Python script

rings_py = py.ring_statistics.ring_statistics(num_atoms_py, coords_py, cutoff, max_ring_size);

% Convert Python result back to MATLAB format
rings = cellfun(@(ring) cellfun(@double, cell(ring)), cell(rings_py), 'UniformOutput', false);

% Display detected rings
disp('Detected Rings:');
for i = 1:length(rings)
    fprintf('Ring %d: %s\n', i, sprintf('%d ', rings{i}));
end

% Plot the ring size distribution
plot_ring_distribution(rings);
saveas(gcf, 'log_files_plot\ring_statistics_plot.png');
% Calculate triatic order parameter q3 for each atom
q3_values = zeros(N, 1);
for i = 1:N
    q3_values(i) = calculate_q3(i, vertex, H1, Lx, Ly);
end

% Calculate the average q3
q3_avg = mean(q3_values);

% Display the average q3
fprintf('Average q3: %f\n', q3_avg);

% Helper function to plot ring size distribution
function plot_ring_distribution(rings)
    ring_sizes = cellfun(@length, rings);
    ring_sizes = ring_sizes(ring_sizes >= 4 & ring_sizes <= 8); % Filter rings with sizes from 4 to 8
    unique_sizes = unique(ring_sizes);
    ring_size_counts = histc(ring_sizes, unique_sizes);
    
    figure;
    bar(unique_sizes, ring_size_counts);
    xlabel('Ring Size');
    ylabel('Number of Rings');
    title('Distribution of Ring Sizes (4 to 8)');
    set(gca, 'XTick', unique_sizes);
    grid on;
end

% Function to calculate triatic order parameter q3 for a given atom
function q3 = calculate_q3(i, vertex, H1, Lx, Ly)
    neighbors = find(H1(i, :));
    num_neighbors = length(neighbors);
    if num_neighbors == 0
        q3 = 0;
        return;
    end

    if num_neighbors == 3
        % Ylm function for l=3 (triatic order)
        l = 3;
        q3_sum = 0; % Initialize sum for q3
        q3=0
     for m = -l:l
        for j = 1:num_neighbors
            r_ij = vertex(neighbors(j), :) - vertex(i, :);
            r_ij = r_ij - round(r_ij ./ [Lx, Ly]) .* [Lx, Ly]; % Periodic boundary conditions
            
            % Convert 2D vector to 3D by adding a zero z-component
            r_ij = [r_ij, 0];
            
            % Compute angles for spherical harmonics
            theta = atan2(norm(cross([1, 0, 0], r_ij)), dot([1, 0, 0], r_ij)); % Angle with respect to x-axis
            phi = atan2(r_ij(2), r_ij(1)); % Azimuthal angle
            P_lm = legendre(l, cos(theta), 'norm');
            Ylm = sqrt((2*l+1)/(4*pi) * factorial(l-abs(m))/factorial(l+abs(m))) * P_lm(abs(m)+1) * exp(1i*m*phi);         
            q3_sum = q3_sum + Ylm;
        end
        q3_sum=q3_sum/num_neighbors;
        q3=q3+abs(q3_sum)^2;

     end
    end

    % Correctly normalize q3
    q3 = sqrt((4 * pi / (2*l + 1)) * q3);
end

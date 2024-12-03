%author: Hassan Shoaib
%email: hassan.shoaib@mail.mcgill.ca

% Create directory for log files
output_dir = 'log_files';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
else
    % If directory exists, delete all files within it
    files = dir(fullfile(output_dir, '*'));
    for i = 1:length(files)
        if ~files(i).isdir
            delete(fullfile(output_dir, files(i).name));
        end
    end
end

% Create directory for log files
output_dir = 'log_files';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Initialize lattice parameters
Nc = 4;   % Number of atoms in a unit cell
Na1 = 70;  % Number of unit cells along x-axis
Na2 = 110; % Number of unit cells along y-axis
Np = Na1 * Na2 * Nc;   % Total number of atoms in the lattice
Nb = ones(Np, 1);   % Initialize neighbor count for each atom
Nb = 3 * Nb;         % Set initial neighbor count to 3 for each atom (honeycomb lattice)
Lx = 3.0 * Na1;      % Lattice width
Ly = sqrt(3) * Na2;  % Lattice height
Nd = 308;     % Number of defects to introduce

% Define the relative positions of atoms within the unit cell
a = [0.5 sqrt(3)/2.0; 1.0 0.0; 2.0 0.0; 2.5 sqrt(3)/2.0];
% Allocate array for storing atom positions
b = zeros(Np, 2);

% Populate the lattice with atoms
for i = 0:Na1-1
    for j = 0:Na2-1
        b((Na1*j+i)*Nc+1:(Na1*j+i)*Nc+Nc, 1) = a(1:Nc, 1) + 3.0 * i;
        b((Na1*j+i)*Nc+1:(Na1*j+i)*Nc+Nc, 2) = a(1:Nc, 2) + sqrt(3) * j;
    end
end

% Apply scaling factor to atom positions
scaling_factor = 1.428; % Adjust this value as needed
b = b * scaling_factor;
Lx = Lx * scaling_factor;
Ly = Ly * scaling_factor;

% Initialize connectivity matrix H to zeros
H = zeros(Np);
% Determine connectivity based on distance between atoms
for i = 1:Np
    for j = 1:i
        dx = abs(b(i, 1) - b(j, 1));
        dy = abs(b(i, 2) - b(j, 2));
        % Apply periodic boundary conditions
        if dx >= Lx/2.0
            dx = Lx - dx;
        end
        if dy >= Ly/2.0
            dy = Ly - dy;
        end
        % If distance close to 1.0, atoms are considered connected
        if abs(sqrt(dx*dx + dy*dy) - scaling_factor) < 1.0e-4
            H(i, j) = 1;
            H(j, i) = 1;
        end
    end
end

% Output the initial connectivity matrix and atom positions
dlmwrite(fullfile(output_dir, 'connectivity_matrix.txt'), H, '\t');

b2(1, 1) = Np;
b2(2, 1) = Lx; b2(2, 2) = 0.0;
b2(3, 1) = 0.0; b2(3, 2) = Ly;
b2(4:3+Np, :) = b(1:Np, :);
dlmwrite(fullfile(output_dir, 'vertex.txt'), b2, 'delimiter', '\t', 'precision', 15);

% Create defects in the lattice
kB = 8.617e-5; % Boltzmann constant
T = 300;    % Temperature
beta = 1/(kB*T);   % Inverse temperature

for n = 1:Nd
    suc = 0;
    while suc == 0
        i = randi(Np, 1);  % Select a random atom
        j2 = randi(Nb(i, 1), 1);    % Randomly select one of its bonds
        t1 = 0;
        j = 1;
        % Find the j-th connected atom
        while t1 < j2 && j <= Np
            if H(i, j) == 1
                t1 = t1 + 1;
            end
            j = j + 1;
        end
        j = j - 1;
        
        % Ensure j is within valid bounds
        if j > Np || j < 1
            continue;
        end

        % Calculate displacement vectors for defect creation
        dx = b(j, 1) - b(i, 1);
        dy = b(j, 2) - b(i, 2);
        % Adjust for periodic boundary
        if dx > Lx/2.0
            dx = dx - Lx;
        elseif dx <= -Lx/2.0
            dx = dx + Lx;
        end
        if dy > Ly/2.0
            dy = dy - Ly;
        elseif dy <= -Ly/2.0
            dy = dy + Ly;
        end
        dx1 = -dy;
        dy1 = dx;
        % Calculate new positions for the defect
        x1 = b(i, 1) + dx1/2.0 + dx/2.0;
        y1 = b(i, 2) + dy1/2.0 + dy/2.0;
        x2 = b(j, 1) - dx1/2.0 - dx/2.0;
        y2 = b(j, 2) - dy1/2.0 - dy/2.0;
        
        % Correct new positions for periodic boundary conditions
        if x1 < 0
            x1 = x1 + Lx;
        elseif x1 >= Lx
            x1 = x1 - Lx;
        end
        if y1 < 0
            y1 = y1 + Ly;
        elseif y1 >= Ly
            y1 = y1 - Ly;
        end
        
        if x2 < 0
            x2 = x2 + Lx;
        elseif x2 >= Lx
            x2 = x2 - Lx;
        end
        if y2 < 0
            y2 = y2 + Ly;
        elseif y2 >= Ly
            y2 = y2 - Ly;
        end
        
        b_temp = b;  % Initialize temporary positions
        b_temp(i, 1) = x1;
        b_temp(i, 2) = y1;
        b_temp(j, 1) = x2;
        b_temp(j, 2) = y2;

        H1 = H; % Copy connectivity matrix for modification

        % Re-evaluate connections to introduce the defect
        % reassigning bonds to create the defect
        for nt = 1:Np
            if (H1(i, nt) == 1) && (nt ~= j)
                dx1 = abs(x1 - b(nt, 1));
                dy1 = abs(y1 - b(nt, 2));
                if dx1 >= Lx / 2.0
                    dx1 = Lx - dx1;
                end
                if dy1 >= Ly / 2.0
                    dy1 = Ly - dy1;
                end
                dist1 = sqrt(dx1 * dx1 + dy1 * dy1);
                
                dx2 = abs(x2 - b(nt, 1));
                dy2 = abs(y2 - b(nt, 2));
                if dx2 >= Lx / 2.0
                    dx2 = Lx - dx2;
                end
                if dy2 >= Ly / 2.0
                    dy2 = Ly - dy2;
                end
                dist2 = sqrt(dx2 * dx2 + dy2 * dy2);
                
                if dist1 > dist2
                    H1(i, nt) = 0;
                    H1(nt, i) = 0;
                    H1(j, nt) = 1;
                    H1(nt, j) = 1;
                end
            end
        end
        
        for nt = 1:Np
            if (H1(j, nt) == 1) && (nt ~= i)
                dx1 = abs(x1 - b(nt, 1));
                dy1 = abs(y1 - b(nt, 2));
                if dx1 >= Lx / 2.0
                    dx1 = Lx - dx1;
                end
                if dy1 >= Ly / 2.0
                    dy1 = Ly - dy1;
                end
                dist1 = sqrt(dx1 * dx1 + dy1 * dy1);
                
                dx2 = abs(x2 - b(nt, 1));
                dy2 = abs(y2 - b(nt, 2));
                if dx2 >= Lx / 2.0
                    dx2 = Lx - dx2;
                end
                if dy2 >= Ly / 2.0
                    dy2 = Ly - dy2;
                end
                dist2 = sqrt(dx2 * dx2 + dy2 * dy2);
                
                if dist1 < dist2
                    H1(i, nt) = 1;
                    H1(nt, i) = 1;
                    H1(j, nt) = 0;
                    H1(nt, j) = 0;
                end
            end
        end

        % Calculate the energies
        E_old = calculate_energy(H, b, Np, Lx, Ly, scaling_factor);
        E_new = calculate_energy(H1, b_temp, Np, Lx, Ly, scaling_factor);
        
        delta_E = E_new - E_old;
        
        % Print debugging information
        fprintf('Energy before: %f, Energy after: %f, ΔE: %f\n', E_old, E_new, delta_E);
        
        % Calculate the acceptance probability
        if delta_E <= 0
            w = 1;
        else
            w = exp((-5*(delta_E-0.06)));
        end
        
        % Generate a random number z
        z = rand;
        
        % Accept the move if w > z
        if w > z && (sum(H1(i, :)) == 3) && (sum(H1(j, :)) == 3)
            b = b_temp;  % Update the positions
            H = H1;     % Update the connectivity matrix
            suc = 1;    % Mark success
            fprintf('Move no. %d Accepted for atoms %d and %d with delta_E = %f\n', n, i, j, delta_E);
        elseif w < z && (sum(H1(i, :)) == 3) && (sum(H1(j, :)) == 3)
            fprintf('Move no. %d rejected for atoms %d and %d with delta_E = %f due to energy\n', n, i, j, delta_E);
        else 
            fprintf('Move no. %d rejected for atoms %d and %d with delta_E = %f due to connectivity\n', n, i, j, delta_E);
        end   
    end
end


% Output the final connectivity matrix and atom positions
dlmwrite(fullfile(output_dir, 'connectivity_matrix.txt'), H, '\t');

b2(1, 1) = Np;
b2(2, 1) = Lx; b2(2, 2) = 0.0;
b2(3, 1) = 0.0; b2(3, 2) = Ly;
b2(4:3+Np, :) = b(1:Np, :);
dlmwrite(fullfile(output_dir, 'vertex.txt'), b2, 'delimiter', '\t', 'precision', 15);

% Relax the structure
relax();

function generate_lammps_data(filename, b, Lx, Ly)
    fid = fopen(filename, 'w');
    fprintf(fid, 'LAMMPS data file\n\n');
    fprintf(fid, '%d atoms\n', size(b, 1));
    fprintf(fid, '1 atom types\n\n');
    fprintf(fid, '0 %f xlo xhi\n', Lx);
    fprintf(fid, '0 %f ylo yhi\n', Ly);
    fprintf(fid, '0 1 zlo zhi\n\n');
    fprintf(fid, 'Atoms\n\n');
    for i = 1:size(b, 1)
        fprintf(fid, '%d 1 %f %f 0.0\n', i, b(i, 1), b(i, 2));
    end
    fclose(fid);
end

function generate_lammps_input(filename, data_filename)
    fid = fopen(filename, 'w');
    fprintf(fid, 'units metal\n');
    fprintf(fid, 'dimension 3\n');
    fprintf(fid, 'boundary p p f\n');
    fprintf(fid, 'atom_style atomic\n');
    fprintf(fid, 'read_data %s\n', fullfile('log_files', data_filename));
    fprintf(fid, 'mass 1 12.01\n');
    fprintf(fid, 'pair_style tersoff\n');
    fprintf(fid, 'pair_coeff * * C.tersoff C\n');
    fprintf(fid, 'neighbor 2.0 bin\n');
    fprintf(fid, 'neigh_modify delay 0 every 1 check yes\n');
    fprintf(fid, 'timestep        0.00025\n');
    fprintf(fid, 'thermo 1\n');
    fprintf(fid, 'thermo_style custom step pe\n');
    fprintf(fid, 'min_style cg\n');  % Conjugate gradient minimization style
    fprintf(fid, 'min_modify dmax 0.1 line quadratic\n');  % Adjust max displacement and line search type
    fprintf(fid, 'minimize 1.0e-10 1.0e-10 10000 10000\n');  % Adjust the tolerance and max steps
    fprintf(fid, 'run 0\n');
    fclose(fid);
end

function energy = extract_energy_from_log(log_filename)
    fid = fopen(log_filename, 'rt');
    energy = NaN;
    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, 'PotEng')
            data = textscan(fid, '%f%f', 'HeaderLines', 0, 'CollectOutput', true);
            energy = data{1}(1, 2); % Get the initial potential energy value
            break;
        end
    end
    fclose(fid);
end

function relax()
    A = dlmread('log_files/vertex.txt');
    N = A(1,1);
    Lx = A(2,1);
    Ly = A(3,2);
    vertex = A(4:3+N,:);
    H1 = dlmread('log_files/connectivity_matrix.txt');
    ct = 3;
    H = zeros(N, ct);
    for i = 1:N
        t = 1;
        for j = 1:N
            if H1(i, j) == 1
                H(i, t) = j;
                t = t + 1;
            end
        end
    end

    num_move = 1000 * N;
    mv_size = 0.1;
    for nt = 1:num_move
        ix = randi(N, 1);
        dE_old = 0;
        dE_new = 0;
        dr = -mv_size / 2.0 + mv_size * rand(1, 2);
        vn = dr + vertex(ix, :);
        if vn(1, 1) > Lx
            vn(1, 1) = vn(1, 1) - Lx;
        elseif vn(1, 1) < 0
            vn(1, 1) = vn(1, 1) + Lx;
        end
        if vn(1, 2) > Ly
            vn(1, 2) = vn(1, 2) - Ly;
        elseif vn(1, 2) < 0
            vn(1, 2) = vn(1, 2) + Ly;
        end
        
        theta_o = zeros(ct, 1);
        theta_n = zeros(ct, 1);
        for j = 1:ct
            % bond length change
            dx = vertex(ix, 1) - vertex(H(ix, j), 1);
            dy = vertex(ix, 2) - vertex(H(ix, j), 2);
            if dx > Lx / 2.0
                dx = dx - Lx;
            elseif dx <= -Lx / 2.0
                dx = dx + Lx;
            end
            if dy > Ly / 2.0
                dy = dy - Ly;
            elseif dy <= -Ly / 2.0
                dy = dy + Ly;
            end
            dE_old = dE_old + (sqrt(dx * dx + dy * dy) - 1.0)^2;
            
            dx = vn(1, 1) - vertex(H(ix, j), 1);
            dy = vn(1, 2) - vertex(H(ix, j), 2);
            if dx > Lx / 2.0
                dx = dx - Lx;
            elseif dx <= -Lx / 2.0
                dx = dx + Lx;
            end
            if dy > Ly / 2.0
                dy = dy - Ly;
            elseif dy <= -Ly / 2.0
                dy = dy + Ly;
            end
            dE_new = dE_new + (sqrt(dx * dx + dy * dy) - 1.0)^2;
            
            % bond angle change for perturbed site
            if j < ct
                x1o = vertex(H(ix, j), :) - vertex(ix, :);
                x2o = vertex(H(ix, j + 1), :) - vertex(ix, :);
                x1n = vertex(H(ix, j), :) - vn;
                x2n = vertex(H(ix, j + 1), :) - vn;
                if x1o(1, 1) > Lx / 2.0
                    x1o(1, 1) = x1o(1, 1) - Lx;
                elseif x1o(1, 1) <= -Lx / 2.0
                    x1o(1, 1) = x1o(1, 1) + Lx;
                end
                if x1o(1, 2) > Ly / 2.0
                    x1o(1, 2) = x1o(1, 2) - Ly;
                elseif x1o(1, 2) <= -Ly / 2.0
                    x1o(1, 2) = x1o(1, 2) + Ly;
                end
                if x1n(1, 1) > Lx / 2.0
                    x1n(1, 1) = x1n(1, 1) - Lx;
                elseif x1n(1, 1) <= -Lx / 2.0
                    x1n(1, 1) = x1n(1, 1) + Lx;
                end
                if x1n(1, 2) > Ly / 2.0
                    x1n(1, 2) = x1n(1, 2) - Ly;
                elseif x1n(1, 2) <= -Ly / 2.0
                    x1n(1, 2) = x1n(1, 2) + Ly;
                end
                if x2o(1, 1) > Lx / 2.0
                    x2o(1, 1) = x2o(1, 1) - Lx;
                elseif x2o(1, 1) <= -Lx / 2.0
                    x2o(1, 1) = x2o(1, 1) + Lx;
                end
                if x2o(1, 2) > Ly / 2.0
                    x2o(1, 2) = x2o(1, 2) - Ly;
                elseif x2o(1, 2) <= -Ly / 2.0
                    x2o(1, 2) = x2o(1, 2) + Ly;
                end
                if x2n(1, 1) > Lx / 2.0
                    x2n(1, 1) = x2n(1, 1) - Lx;
                elseif x2n(1, 1) <= -Lx / 2.0
                    x2n(1, 1) = x2n(1, 1) + Lx;
                end
                if x2n(1, 2) > Ly / 2.0
                    x2n(1, 2) = x2n(1, 2) - Ly;
                elseif x2n(1, 2) <= -Ly / 2.0
                    x2n(1, 2) = x2n(1, 2) + Ly;
                end
                CosTheta_o = max(min(dot(x1o, x2o) / (norm(x1o) * norm(x2o)), 1), -1);
                CosTheta_n = max(min(dot(x1n, x2n) / (norm(x1n) * norm(x2n)), 1), -1);
                theta_o(j, 1) = acos(CosTheta_o);
                theta_n(j, 1) = acos(CosTheta_n);
            else
                x1o = vertex(H(ix, j), :) - vertex(ix, :);
                x2o = vertex(H(ix, 1), :) - vertex(ix, :);
                x1n = vertex(H(ix, j), :) - vn;
                x2n = vertex(H(ix, 1), :) - vn;
                if x1o(1, 1) > Lx / 2.0
                    x1o(1, 1) = x1o(1, 1) - Lx;
                elseif x1o(1, 1) <= -Lx / 2.0
                    x1o(1, 1) = x1o(1, 1) + Lx;
                end
                if x1o(1, 2) > Ly / 2.0
                    x1o(1, 2) = x1o(1, 2) - Ly;
                elseif x1o(1, 2) <= -Ly / 2.0
                    x1o(1, 2) = x1o(1, 2) + Ly;
                end
                if x1n(1, 1) > Lx / 2.0
                    x1n(1, 1) = x1n(1, 1) - Lx;
                elseif x1n(1, 1) <= -Lx / 2.0
                    x1n(1, 1) = x1n(1, 1) + Lx;
                end
                if x1n(1, 2) > Ly / 2.0
                    x1n(1, 2) = x1n(1, 2) - Ly;
                elseif x1n(1, 2) <= -Ly / 2.0
                    x1n(1, 2) = x1n(1, 2) + Ly;
                end
                
                if x2o(1, 1) > Lx / 2.0
                    x2o(1, 1) = x2o(1, 1) - Lx;
                elseif x2o(1, 1) <= -Lx / 2.0
                    x2o(1, 1) = x2o(1, 1) + Lx;
                end
                if x2o(1, 2) > Ly / 2.0
                    x2o(1, 2) = x2o(1, 2) - Ly;
                elseif x2o(1, 2) <= -Ly / 2.0
                    x2o(1, 2) = x2o(1, 2) + Ly;
                end
                if x2n(1, 1) > Lx / 2.0
                    x2n(1, 1) = x2n(1, 1) - Lx;
                elseif x2n(1, 1) <= -Lx / 2.0
                    x2n(1, 1) = x2n(1, 1) + Lx;
                end
                if x2n(1, 2) > Ly / 2.0
                    x2n(1, 2) = x2n(1, 2) - Ly;
                elseif x2n(1, 2) <= -Ly / 2.0
                    x2n(1, 2) = x2n(1, 2) + Ly;
                end
                CosTheta_o = max(min(dot(x1o, x2o) / (norm(x1o) * norm(x2o)), 1), -1);
                CosTheta_n = max(min(dot(x1n, x2n) / (norm(x1n) * norm(x2n)), 1), -1);
                theta_o(j, 1) = acos(CosTheta_o);
                theta_n(j, 1) = acos(CosTheta_n);
            end
            
        end
        
        if abs(theta_o(1, 1) + theta_o(2, 1) - theta_o(3, 1)) < 1e-8
            theta_o(3, 1) = 2 * pi - theta_o(3, 1);
        end
        if abs(theta_o(1, 1) + theta_o(3, 1) - theta_o(2, 1)) < 1e-8
            theta_o(2, 1) = 2 * pi - theta_o(2, 1);
        end
        if abs(theta_o(2, 1) + theta_o(3, 1) - theta_o(1, 1)) < 1e-8
            theta_o(1, 1) = 2 * pi - theta_o(1, 1);
        end
        
        if abs(theta_n(1, 1) + theta_n(2, 1) - theta_n(3, 1)) < 1e-8
            theta_n(3, 1) = 2 * pi - theta_n(3, 1);
        end
        if abs(theta_n(1, 1) + theta_n(3, 1) - theta_n(2, 1)) < 1e-8
            theta_n(2, 1) = 2 * pi - theta_n(2, 1);
        end
        if abs(theta_n(2, 1) + theta_n(3, 1) - theta_n(1, 1)) < 1e-8
            theta_n(1, 1) = 2 * pi - theta_n(1, 1);
        end
        
        for j = 1:ct
            dE_old = dE_old + (theta_o(j, 1) - pi * 2.0 / 3.0)^2;
            dE_new = dE_new + (theta_n(j, 1) - pi * 2.0 / 3.0)^2;
        end
        
        % bond angle change for neighboring sites
        for j = 1:ct
            jx = H(ix, j);
            theta_o = zeros(ct, 1);
            theta_n = zeros(ct, 1);
            for k = 1:ct
                if k < ct
                    x1o = vertex(H(jx, k), :) - vertex(jx, :);
                    x2o = vertex(H(jx, k + 1), :) - vertex(jx, :);
                    x1n = x1o;
                    x2n = x2o;
                    if H(jx, k) == ix
                        x1n = x1n + dr;
                    end
                    if H(jx, k + 1) == ix
                        x2n = x2n + dr;
                    end
                    if x1o(1, 1) > Lx / 2.0
                        x1o(1, 1) = x1o(1, 1) - Lx;
                    elseif x1o(1, 1) <= -Lx / 2.0
                        x1o(1, 1) = x1o(1, 1) + Lx;
                    end
                    if x1o(1, 2) > Ly / 2.0
                        x1o(1, 2) = x1o(1, 2) - Ly;
                    elseif x1o(1, 2) <= -Ly / 2.0
                        x1o(1, 2) = x1o(1, 2) + Ly;
                    end
                    if x1n(1, 1) > Lx / 2.0
                        x1n(1, 1) = x1n(1, 1) - Lx;
                    elseif x1n(1, 1) <= -Lx / 2.0
                        x1n(1, 1) = x1n(1, 1) + Lx;
                    end
                    if x1n(1, 2) > Ly / 2.0
                        x1n(1, 2) = x1n(1, 2) - Ly;
                    elseif x1n(1, 2) <= -Ly / 2.0
                        x1n(1, 2) = x1n(1, 2) + Ly;
                    end
                    if x2o(1, 1) > Lx / 2.0
                        x2o(1, 1) = x2o(1, 1) - Lx;
                    elseif x2o(1, 1) <= -Lx / 2.0
                        x2o(1, 1) = x2o(1, 1) + Lx;
                    end
                    if x2o(1, 2) > Ly / 2.0
                        x2o(1, 2) = x2o(1, 2) - Ly;
                    elseif x2o(1, 2) <= -Ly / 2.0
                        x2o(1, 2) = x2o(1, 2) + Ly;
                    end
                    if x2n(1, 1) > Lx / 2.0
                        x2n(1, 1) = x2n(1, 1) - Lx;
                    elseif x2n(1, 1) <= -Lx / 2.0
                        x2n(1, 1) = x2n(1, 1) + Lx;
                    end
                    if x2n(1, 2) > Ly / 2.0
                        x2n(1, 2) = x2n(1, 2) - Ly;
                    elseif x2n(1, 2) <= -Ly / 2.0
                        x2n(1, 2) = x2n(1, 2) + Ly;
                    end
                    
                    CosTheta_o = max(min(dot(x1o, x2o) / (norm(x1o) * norm(x2o)), 1), -1);
                    CosTheta_n = max(min(dot(x1n, x2n) / (norm(x1n) * norm(x2n)), 1), -1);
                    theta_o(k, 1) = acos(CosTheta_o);
                    theta_n(k, 1) = acos(CosTheta_n);
                else
                    x1o = vertex(H(jx, k), :) - vertex(jx, :);
                    x2o = vertex(H(jx, 1), :) - vertex(jx, :);
                    x1n = x1o;
                    x2n = x2o;
                    if H(jx, k) == ix
                        x1n = x1n + dr;
                    end
                    if H(jx, 1) == ix
                        x2n = x2n + dr;
                    end
                    if x1o(1, 1) > Lx / 2.0
                        x1o(1, 1) = x1o(1, 1) - Lx;
                    elseif x1o(1, 1) <= -Lx / 2.0
                        x1o(1, 1) = x1o(1, 1) + Lx;
                    end
                    if x1o(1, 2) > Ly / 2.0
                        x1o(1, 2) = x1o(1, 2) - Ly;
                    elseif x1o(1, 2) <= -Ly / 2.0
                        x1o(1, 2) = x1o(1, 2) + Ly;
                    end
                    if x1n(1, 1) > Lx / 2.0
                        x1n(1, 1) = x1n(1, 1) - Lx;
                    elseif x1n(1, 1) <= -Lx / 2.0
                        x1n(1, 1) = x1n(1, 1) + Lx;
                    end
                    if x1n(1, 2) > Ly / 2.0
                        x1n(1, 2) = x1n(1, 2) - Ly;
                    elseif x1n(1, 2) <= -Ly / 2.0
                        x1n(1, 2) = x1n(1, 2) + Ly;
                    end
                    if x2o(1, 1) > Lx / 2.0
                        x2o(1, 1) = x2o(1, 1) - Lx;
                    elseif x2o(1, 1) <= -Lx / 2.0
                        x2o(1, 1) = x2o(1, 1) + Lx;
                    end
                    if x2o(1, 2) > Ly / 2.0
                        x2o(1, 2) = x2o(1, 2) - Ly;
                    elseif x2o(1, 2) <= -Ly / 2.0
                        x2o(1, 2) = x2o(1, 2) + Ly;
                    end
                    if x2n(1, 1) > Lx / 2.0
                        x2n(1, 1) = x2n(1, 1) - Lx;
                    elseif x2n(1, 1) <= -Lx / 2.0
                        x2n(1, 1) = x2n(1, 1) + Lx;
                    end
                    if x2n(1, 2) > Ly / 2.0
                        x2n(1, 2) = x2n(1, 2) - Ly;
                    elseif x2n(1, 2) <= -Ly / 2.0
                        x2n(1, 2) = x2n(1, 2) + Ly;
                    end
                    
                    CosTheta_o = max(min(dot(x1o, x2o) / (norm(x1o) * norm(x2o)), 1), -1);
                    CosTheta_n = max(min(dot(x1n, x2n) / (norm(x1n) * norm(x2n)), 1), -1);
                    theta_o(k, 1) = acos(CosTheta_o);
                    theta_n(k, 1) = acos(CosTheta_n);
                end
            end
            
            if abs(theta_o(1, 1) + theta_o(2, 1) - theta_o(3, 1)) < 1e-8
                theta_o(3, 1) = 2 * pi - theta_o(3, 1);
            end
            if abs(theta_o(1, 1) + theta_o(3, 1) - theta_o(2, 1)) < 1e-8
                theta_o(2, 1) = 2 * pi - theta_o(2, 1);
            end
            if abs(theta_o(2, 1) + theta_o(3, 1) - theta_o(1, 1)) < 1e-8
                theta_o(1, 1) = 2 * pi - theta_o(1, 1);
            end
            
            if abs(theta_n(1, 1) + theta_n(2, 1) - theta_n(3, 1)) < 1e-8
                theta_n(3, 1) = 2 * pi - theta_n(3, 1);
            end
            if abs(theta_n(1, 1) + theta_n(3, 1) - theta_n(2, 1)) < 1e-8
                theta_n(2, 1) = 2 * pi - theta_n(2, 1);
            end
            if abs(theta_n(2, 1) + theta_n(3, 1) - theta_n(1, 1)) < 1e-8
                theta_n(1, 1) = 2 * pi - theta_n(1, 1);
            end
            
            for k = 1:ct
                dE_old = dE_old + (theta_o(k, 1) - pi * 2.0 / 3.0)^2;
                dE_new = dE_new + (theta_n(k, 1) - pi * 2.0 / 3.0)^2;
            end
        end
        
        if dE_old > dE_new
            vertex(ix, :) = vn;
        end
    end

    A(4:3+N, :) = vertex;
    dlmwrite(fullfile('log_files', 'vertex.txt'), A, 'delimiter', '\t', 'precision', 15);
end

function b = read_lammps_data(filename)
    fid = fopen(filename, 'r');
    tline = fgetl(fid);
    while ischar(tline)
        if contains(tline, 'Atoms')
            break;
        end
        tline = fgetl(fid);
    end
    tline = fgetl(fid); % Skip the 'Atoms' line
    data = textscan(fid, '%d %d %f %f %f');
    b = [data{3} data{4}];
    fclose(fid);
end

function H = update_connectivity_matrix(b, Np, Lx, Ly, scaling_factor)
    H = zeros(Np);
    for i = 1:Np
        for j = 1:i
            dx = abs(b(i, 1) - b(j, 1));
            dy = abs(b(i, 2) - b(j, 2));
            if dx >= Lx / 2.0
                dx = Lx - dx;
            end
            if dy >= Ly / 2.0
                dy = Ly - dy;
            end
            if abs(sqrt(dx * dx + dy * dy) - scaling_factor) < 1.0e-4
                H(i, j) = 1;
                H(j, i) = 1;
            end
        end
    end
end

function energy = calculate_energy(H, b, Np, Lx, Ly, scaling_factor)
    energy = 0;
    for i = 1:Np
        for j = 1:Np
            if H(i, j) == 1
                dx = abs(b(i, 1) - b(j, 1));
                dy = abs(b(i, 2) - b(j, 2));
                % Apply periodic boundary conditions
                if dx > Lx / 2.0
                    dx = dx - Lx;
                elseif dx <= -Lx / 2.0
                    dx = dx + Lx;
                end
                if dy > Ly / 2.0
                    dy = dy - Ly;
                elseif dy <= -Ly / 2.0
                    dy = dy + Ly;
                end
                dist = sqrt(dx^2 + dy^2);
                energy = energy + (dist - scaling_factor)^2;
            end
        end
    end
end

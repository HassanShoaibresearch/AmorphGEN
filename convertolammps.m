% Load the vertex data
vertexData = dlmread('log_files/vertex.txt', '\t');
N = vertexData(1, 1); % Assuming the number of atoms is the first entry
positions = vertexData(4:end, :); % Assuming atom positions start at line 4

% Scaling factor to adjust the distances to match graphene's bond length
scalingFactor = 1;

% Apply the scaling factor to the atom positions
scaledPositions = positions * scalingFactor;

% Optionally, you might need to adjust Lx and Ly if they are part of the output
Lx = vertexData(2, 1) * scalingFactor;
Ly = vertexData(3, 2) * scalingFactor;

scaledvertex(1,1)=N;
scaledvertex(2,1)=Lx;scaledvertex(2,2)=0.0;
scaledvertex(3,1)=0.0;scaledvertex(3,2)=Ly;
scaledvertex(4:3+N,:)=scaledPositions(1:N,:);

% Write the scaled vertex data back to a new file
dlmwrite('./vertex.txt', scaledvertex, 'delimiter', '\t', 'precision', 15);

% Load vertex data
vertexData = dlmread('./vertex.txt', '\t');
N = vertexData(1, 1); % Number of atoms
Lx = vertexData(2, 1);
Ly = vertexData(3, 2);
atomPositions = vertexData(4:end, :);

% Assuming one atom type for all
atomType = 1;

% Assuming connectivity matrix is available and named 'connectivity_matrix.txt'
% This part is optional and depends on whether you want to include bonds
connectivity = dlmread('./connectivity_matrix.txt', '\t');
[bondList, ~] = find(tril(connectivity)); % Extracting lower triangular to avoid duplicates

% Preparing to write to LAMMPS data file
dataFileName = './lammps_data_file.data';
fid = fopen(dataFileName, 'w');

% Header
fprintf(fid, '# LAMMPS data file generated from vertex.txt\n\n');
fprintf(fid, '%d atoms\n', N);
fprintf(fid, '%d bonds\n\n', length(bondList)/2);

fprintf(fid, '1 atom types\n');
if ~isempty(bondList)
    fprintf(fid, '1 bond types\n\n');
end

% Box dimensions
fprintf(fid, '0.0 %f xlo xhi\n', Lx);
fprintf(fid, '0.0 %f ylo yhi\n', Ly);
fprintf(fid, '-100 100 zlo zhi\n\n');

% Masses section (assuming a mass of 1.0 for simplicity)
fprintf(fid, 'Masses\n\n');
fprintf(fid, '1 1.0\n\n');

% Atoms section
fprintf(fid, 'Atoms\n\n');
for i = 1:N
    fprintf(fid, '%d %d %f %f %f\n', i, atomType, atomPositions(i, 1), atomPositions(i, 2), 0.0);
end

% Bonds section
if ~isempty(bondList)
    fprintf(fid, '\nBonds\n\n');
    for i = 1:2:length(bondList) % Skipping every other to avoid duplicates
        fprintf(fid, '%d 1 %d %d\n', (i+1)/2, bondList(i), bondList(i+1));
    end
end

fclose(fid);

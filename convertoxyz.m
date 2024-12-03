% Load data from vertex file
fileID = fopen('log_files\vertex.txt', 'r');
data = textscan(fileID, '%f %f');
fclose(fileID);

% Extract number of atoms, x-axis length, and y-axis length
num_atoms = data{1}(1);
Lx = data{1}(2);
Ly = data{2}(3);

% Extract atom positions
positions = [data{1}(4:end), data{2}(4:end)];

% Open the new XYZ file for writing
fileID = fopen('log_files_plot\model.xyz', 'w');

% Write the number of atoms
fprintf(fileID, '%d\n', num_atoms);

% Write the header line with lattice parameters
fprintf(fileID, 'pbc="T T F" Lattice="%f 0 0 0 %f 0 0 0 3.35" Properties=species:S:1:pos:R:3:group:I:1\n', Lx, Ly);

% Write the atom positions
for i = 1:num_atoms
    fprintf(fileID, 'C %f %f 0 %d\n', positions(i, 1), positions(i, 2),i-1);
end

% Close the XYZ file
fclose(fileID);
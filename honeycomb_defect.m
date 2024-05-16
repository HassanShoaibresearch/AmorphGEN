%author: Hassan Shoaib
%email: hassan.shoaib@mail.mcgill.ca

% Initialize lattice parameters
Nc=4;   % Number of atoms in a unit cell
Na1=35;  % Number of unit cells along x-axis
Na2=55; % Number of unit cells along y-axis
Np=Na1*Na2*Nc   % Total number of atoms in the lattice
Nb=ones(Np,1);  % Initialize neighbor count for each atom
Nb=3*Nb;        % Set initial neighbor count to 3 for each atom (honeycomb lattice)
Lx=3.0*Na1;     % Lattice width
Ly=sqrt(3)*Na2; % Lattice height


% Define the relative positions of atoms within the unit cell
a=[0.5 sqrt(3)/2.0; 1.0 0.0; 2.0 0.0; 2.5 sqrt(3)/2.0];
% Allocate array for storing atom positions
b=zeros(Np,2);

% Populate the lattice with atoms
for i=0:Na1-1
    for j=0:Na2-1
        b((Na1*j+i)*Nc+1:(Na1*j+i)*Nc+Nc,1)=a(1:Nc,1)+3.0*i;
        b((Na1*j+i)*Nc+1:(Na1*j+i)*Nc+Nc,2)=a(1:Nc,2)+sqrt(3)*j;
    end
end

% Initialize connectivity matrix H to zeros
H=zeros(Np);
% Determine connectivity based on distance between atoms
for i=1:Np
    for j=1:i
        dx=abs(b(i,1)-b(j,1));
        dy=abs(b(i,2)-b(j,2));
        % Apply periodic boundary conditions
        if dx >= Lx/2.0
            dx = Lx-dx;
        end
        if dy >= Ly/2.0
            dy = Ly-dy;
        end
        % If distance close to 1.0, atoms are considered connected
        if abs(sqrt(dx*dx+dy*dy)-1.0)<1.0e-4
            H(i,j)=1;H(j,i)=1;
        end
    end
end

% Create defects in the lattice
Nd = 30;     % Number of defects to introduce

%mov_ind=zeros(Np,1);
%b1=zeros(Nd+Np,2);
%H1=zeros(Nd+Np);
%b1(1:Np,1:2)=b(1:Np,1:2);
%H1(1:Np,1:Np)=H(1:Np,1:Np);
%theta=0.4;

for n=1:Nd
    suc = 0;
    while suc == 0
    i=randi(Np,1);  % Select a random atom
    %while mov_ind(i,1)==1
        %i=randi(Np,1);
    %end
    j2=randi(Nb(i,1),1);    % Randomly select one of its bonds
    t1=0;
    j=1;
    % Find the j-th connected atom
    while t1<j2
        if H(i,j)==1
            t1=t1+1;
        end
        j=j+1;
    end
    j=j-1;
    % Calculate displacement vectors for defect creation
    dx=b(j,1)-b(i,1);
    dy=b(j,2)-b(i,2);
    % Adjust for periodic boundary
    if dx > Lx/2.0
        dx=dx-Lx;
    elseif dx <= -Lx/2.0
        dx=dx+Lx;
    end
    if dy > Ly/2.0
        dy=dy-Ly;
    elseif dy <= -Ly/2.0
        dy=dy+Ly;
    end
    dx1=-dy;
    dy1=dx;
    % Calculate new positions for the defect
    x1=b(i,1)+dx1/2.0+dx/2.0;
    y1=b(i,2)+dy1/2.0+dy/2.0;
    x2=b(j,1)-dx1/2.0-dx/2.0;
    y2=b(j,2)-dy1/2.0-dy/2.0;
    
    % Correct new positions for periodic boundary conditions
    if x1<0
        x1=x1+Lx;
    elseif x1>=Lx
        x1=x1-Lx;
    end
    if y1<0
        y1=y1+Ly;
    elseif y1>=Ly
        y1=y1-Ly;
    end
    
    if x2<0
        x2=x2+Lx;
    elseif x2>=Lx
        x2=x2-Lx;
    end
    if y2<0
        y2=y2+Ly;
    elseif y2>=Ly
        y2=y2-Ly;
    end
    
    H1 = H; % Copy connectivity matrix for modification

    % Re-evaluate connections to introduce the defect
    % reassigning bonds to create the defect
    for nt = 1:Np
        if (H1(i,nt) == 1) && (nt~=j)
            dx1=abs(x1-b(nt,1));
            dy1=abs(y1-b(nt,2));
            if dx1 >= Lx/2.0
                dx1 = Lx-dx1;
            end
            if dy1 >= Ly/2.0
                dy1 = Ly-dy1;
            end
            dist1 = sqrt(dx1*dx1+dy1*dy1);
            
            dx2=abs(x2-b(nt,1));
            dy2=abs(y2-b(nt,2));
            if dx2 >= Lx/2.0
                dx2 = Lx-dx2;
            end
            if dy2 >= Ly/2.0
                dy2 = Ly-dy2;
            end
            dist2 = sqrt(dx2*dx2+dy2*dy2);
            
            if dist1 > dist2
                H1(i,nt) = 0;
                H1(nt,i) = 0;
                H1(j,nt) = 1;
                H1(nt,j) = 1;
            end
        end
    end
    
    for nt = 1:Np
        if (H1(j,nt) == 1) && (nt~=i)
            dx1=abs(x1-b(nt,1));
            dy1=abs(y1-b(nt,2));
            if dx1 >= Lx/2.0
                dx1 = Lx-dx1;
            end
            if dy1 >= Ly/2.0
                dy1 = Ly-dy1;
            end
            dist1 = sqrt(dx1*dx1+dy1*dy1);
            
            dx2=abs(x2-b(nt,1));
            dy2=abs(y2-b(nt,2));
            if dx2 >= Lx/2.0
                dx2 = Lx-dx2;
            end
            if dy2 >= Ly/2.0
                dy2 = Ly-dy2;
            end
            dist2 = sqrt(dx2*dx2+dy2*dy2);
            
            if dist1 < dist2
                H1(i,nt) = 1;
                H1(nt,i) = 1;
                H1(j,nt) = 0;
                H1(nt,j) = 0;
            end
        end
    end

    % Check if defect creation was successful and update positions and connectivity
    if (sum(H1(i,:))==3) && (sum(H1(j,:))==3)
        b(i,1) = x1;
        b(i,2) = y1;
        b(j,1) = x2;
        b(j,2) = y2;
        H = H1;     % Update the connectivity matrix
        suc = 1;    % Mark success
    end
    end
end

%{
for i = 1:Np
    if sum(H(i,:)) > 3
        i
        k = 0;
        for j = 1:Np
            if (H(i,j) == 1) && (sum(H(j,:)) < 3)
                k = j;
            end
        end
        jt = 0;
        min_dist = 1000;
        for j = 1:Np
            if (j~=k) && (H(i,j)==1)
                dx2=abs(b(j,1)-b(k,1));
                dy2=abs(b(j,2)-b(k,2));
                if dx2 >= Lx/2.0
                    dx2 = Lx-dx2;
                end
                if dy2 >= Ly/2.0
                    dy2 = Ly-dy2;
                end
                dist2 = sqrt(dx2*dx2+dy2*dy2);
                if dist2 < min_dist
                    jt = j;
                    min_dist = dist2;
                end
            end
        end
        H(i,jt) = 0;
        H(jt,i) = 0;
        H(k,jt) = 1;
        H(jt,k) = 1;
    end
end
%}

% Output the final connectivity matrix and atom positions
dlmwrite('./connectivity_matrix.txt', H, '\t');

b2(1,1)=Np;
b2(2,1)=Lx;b2(2,2)=0.0;
b2(3,1)=0.0;b2(3,2)=Ly;
b2(4:3+Np,:)=b(1:Np,:);
dlmwrite('./vertex.txt', b2, 'delimiter','\t','precision',15);
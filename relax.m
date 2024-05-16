A=dlmread('./vertex.txt');
%A=dlmread('./vertex_honeycomb.txt');
N=A(1,1);
Lx=A(2,1);
Ly=A(3,2);
vertex=A(4:3+N,:);
H1=dlmread('./connectivity_matrix.txt');
%H1=dlmread('./connectivity_matrix_honeycomb.txt');
ct = 3;
H=zeros(N,ct);
for i = 1:N
    t = 1;
    for j = 1:N
        if H1(i,j) == 1
            H(i,t) = j;
            t = t + 1;
        end
    end
end

num_move = 1000*N;
mv_size = 0.1;
dE = zeros(num_move,1);
for nt = 1:num_move
    ix = randi(N, 1);
    dE_old = 0;
    dE_new = 0;
    dr = -mv_size / 2.0 + mv_size * rand(1,2);
    vn = dr + vertex(ix,:);
    if vn(1,1) > Lx
        vn(1,1) = vn(1,1) - Lx;
    elseif vn(1,1) < 0
        vn(1,1) = vn(1,1) + Lx;
    end
    if vn(1,2) > Ly
        vn(1,2) = vn(1,2) - Ly;
    elseif vn(1,2) < 0
        vn(1,2) = vn(1,2) + Ly;
    end
    
    theta_o = zeros(ct,1);
    theta_n = zeros(ct,1);
    for j = 1:ct
        %bond length change
        dx = vertex(ix,1) - vertex(H(ix,j),1);
        dy = vertex(ix,2) - vertex(H(ix,j),2);
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
        dE_old = dE_old + (sqrt(dx*dx+dy*dy) - 1.0)^2;
        
        dx = vn(1,1) - vertex(H(ix,j),1);
        dy = vn(1,2) - vertex(H(ix,j),2);
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
        dE_new = dE_new + (sqrt(dx*dx+dy*dy) - 1.0)^2;
        
        %bond angle change for perturbed site
        if j<ct
            x1o = vertex(H(ix,j),:) - vertex(ix,:);
            x2o = vertex(H(ix,j+1),:) - vertex(ix,:);
            x1n = vertex(H(ix,j),:) - vn;
            x2n = vertex(H(ix,j+1),:) - vn;
            if x1o(1,1) > Lx/2.0
                x1o(1,1)=x1o(1,1)-Lx;
            elseif x1o(1,1) <= -Lx/2.0
                x1o(1,1)=x1o(1,1)+Lx;
            end
            if x1o(1,2) > Ly/2.0
                x1o(1,2)=x1o(1,2)-Ly;
            elseif x1o(1,2) <= -Ly/2.0
                x1o(1,2)=x1o(1,2)+Ly;
            end
            if x1n(1,1) > Lx/2.0
                x1n(1,1)=x1n(1,1)-Lx;
            elseif x1n(1,1) <= -Lx/2.0
                x1n(1,1)=x1n(1,1)+Lx;
            end
            if x1n(1,2) > Ly/2.0
                x1n(1,2)=x1n(1,2)-Ly;
            elseif x1n(1,2) <= -Ly/2.0
                x1n(1,2)=x1n(1,2)+Ly;
            end
            if x2o(1,1) > Lx/2.0
                x2o(1,1)=x2o(1,1)-Lx;
            elseif x2o(1,1) <= -Lx/2.0
                x2o(1,1)=x2o(1,1)+Lx;
            end
            if x2o(1,2) > Ly/2.0
                x2o(1,2)=x2o(1,2)-Ly;
            elseif x2o(1,2) <= -Ly/2.0
                x2o(1,2)=x2o(1,2)+Ly;
            end
            if x2n(1,1) > Lx/2.0
                x2n(1,1)=x2n(1,1)-Lx;
            elseif x2n(1,1) <= -Lx/2.0
                x2n(1,1)=x2n(1,1)+Lx;
            end
            if x2n(1,2) > Ly/2.0
                x2n(1,2)=x2n(1,2)-Ly;
            elseif x2n(1,2) <= -Ly/2.0
                x2n(1,2)=x2n(1,2)+Ly;
            end
            CosTheta_o = max(min(dot(x1o,x2o)/(norm(x1o)*norm(x2o)),1),-1);
            CosTheta_n = max(min(dot(x1n,x2n)/(norm(x1n)*norm(x2n)),1),-1);
            theta_o(j,1) = acos(CosTheta_o);
            theta_n(j,1) = acos(CosTheta_n);
        else
            x1o = vertex(H(ix,j),:) - vertex(ix,:);
            x2o = vertex(H(ix,1),:) - vertex(ix,:);
            x1n = vertex(H(ix,j),:) - vn;
            x2n = vertex(H(ix,1),:) - vn;
            if x1o(1,1) > Lx/2.0
                x1o(1,1)=x1o(1,1)-Lx;
            elseif x1o(1,1) <= -Lx/2.0
                x1o(1,1)=x1o(1,1)+Lx;
            end
            if x1o(1,2) > Ly/2.0
                x1o(1,2)=x1o(1,2)-Ly;
            elseif x1o(1,2) <= -Ly/2.0
                x1o(1,2)=x1o(1,2)+Ly;
            end
            if x1n(1,1) > Lx/2.0
                x1n(1,1)=x1n(1,1)-Lx;
            elseif x1n(1,1) <= -Lx/2.0
                x1n(1,1)=x1n(1,1)+Lx;
            end
            if x1n(1,2) > Ly/2.0
                x1n(1,2)=x1n(1,2)-Ly;
            elseif x1n(1,2) <= -Ly/2.0
                x1n(1,2)=x1n(1,2)+Ly;
            end
            
            if x2o(1,1) > Lx/2.0
                x2o(1,1)=x2o(1,1)-Lx;
            elseif x2o(1,1) <= -Lx/2.0
                x2o(1,1)=x2o(1,1)+Lx;
            end
            if x2o(1,2) > Ly/2.0
                x2o(1,2)=x2o(1,2)-Ly;
            elseif x2o(1,2) <= -Ly/2.0
                x2o(1,2)=x2o(1,2)+Ly;
            end
            if x2n(1,1) > Lx/2.0
                x2n(1,1)=x2n(1,1)-Lx;
            elseif x2n(1,1) <= -Lx/2.0
                x2n(1,1)=x2n(1,1)+Lx;
            end
            if x2n(1,2) > Ly/2.0
                x2n(1,2)=x2n(1,2)-Ly;
            elseif x2n(1,2) <= -Ly/2.0
                x2n(1,2)=x2n(1,2)+Ly;
            end
            CosTheta_o = max(min(dot(x1o,x2o)/(norm(x1o)*norm(x2o)),1),-1);
            CosTheta_n = max(min(dot(x1n,x2n)/(norm(x1n)*norm(x2n)),1),-1);
            theta_o(j,1) = acos(CosTheta_o);
            theta_n(j,1) = acos(CosTheta_n);
        end
        
    end
    
     if abs(theta_o(1,1) + theta_o(2,1) - theta_o(3,1)) < 1e-8
         theta_o(3,1) = 2*pi - theta_o(3,1);
     end
     if abs(theta_o(1,1) + theta_o(3,1) - theta_o(2,1)) < 1e-8
         theta_o(2,1) = 2*pi - theta_o(2,1);
     end
     if abs(theta_o(2,1) + theta_o(3,1) - theta_o(1,1)) < 1e-8
         theta_o(1,1) = 2*pi - theta_o(1,1);
     end
    
     if abs(theta_n(1,1) + theta_n(2,1) - theta_n(3,1)) < 1e-8
         theta_n(3,1) = 2*pi - theta_n(3,1);
     end
     if abs(theta_n(1,1) + theta_n(3,1) - theta_n(2,1)) < 1e-8
         theta_n(2,1) = 2*pi - theta_n(2,1);
     end
     if abs(theta_n(2,1) + theta_n(3,1) - theta_n(1,1)) < 1e-8
         theta_n(1,1) = 2*pi - theta_n(1,1);
     end
    
    for j = 1:ct
        dE_old = dE_old + (theta_o(j,1) - pi*2.0/3.0)^2;
        dE_new = dE_new + (theta_n(j,1) - pi*2.0/3.0)^2;
        
    end
    
    %bond angle change for neighboring sites
    for j = 1:ct
        jx = H(ix,j);
        theta_o = zeros(ct,1);
        theta_n = zeros(ct,1);
        for k = 1:ct
            if k<ct
                x1o = vertex(H(jx,k),:) - vertex(jx,:);
                x2o = vertex(H(jx,k+1),:) - vertex(jx,:);
                x1n = x1o;
                x2n = x2o;
                if H(jx,k) == ix
                    x1n = x1n + dr;
                end
                if H(jx,k+1) == ix
                    x2n = x2n + dr;
                end
                if x1o(1,1) > Lx/2.0
                    x1o(1,1)=x1o(1,1)-Lx;
                elseif x1o(1,1) <= -Lx/2.0
                    x1o(1,1)=x1o(1,1)+Lx;
                end
                if x1o(1,2) > Ly/2.0
                    x1o(1,2)=x1o(1,2)-Ly;
                elseif x1o(1,2) <= -Ly/2.0
                    x1o(1,2)=x1o(1,2)+Ly;
                end
                if x1n(1,1) > Lx/2.0
                    x1n(1,1)=x1n(1,1)-Lx;
                elseif x1n(1,1) <= -Lx/2.0
                    x1n(1,1)=x1n(1,1)+Lx;
                end
                if x1n(1,2) > Ly/2.0
                    x1n(1,2)=x1n(1,2)-Ly;
                elseif x1n(1,2) <= -Ly/2.0
                    x1n(1,2)=x1n(1,2)+Ly;
                end
                if x2o(1,1) > Lx/2.0
                    x2o(1,1)=x2o(1,1)-Lx;
                elseif x2o(1,1) <= -Lx/2.0
                    x2o(1,1)=x2o(1,1)+Lx;
                end
                if x2o(1,2) > Ly/2.0
                    x2o(1,2)=x2o(1,2)-Ly;
                elseif x2o(1,2) <= -Ly/2.0
                    x2o(1,2)=x2o(1,2)+Ly;
                end
                if x2n(1,1) > Lx/2.0
                    x2n(1,1)=x2n(1,1)-Lx;
                elseif x2n(1,1) <= -Lx/2.0
                    x2n(1,1)=x2n(1,1)+Lx;
                end
                if x2n(1,2) > Ly/2.0
                    x2n(1,2)=x2n(1,2)-Ly;
                elseif x2n(1,2) <= -Ly/2.0
                    x2n(1,2)=x2n(1,2)+Ly;
                end
                
                CosTheta_o = max(min(dot(x1o,x2o)/(norm(x1o)*norm(x2o)),1),-1);
                CosTheta_n = max(min(dot(x1n,x2n)/(norm(x1n)*norm(x2n)),1),-1);
                theta_o(k,1) = acos(CosTheta_o);
                theta_n(k,1) = acos(CosTheta_n);
            else
                x1o = vertex(H(jx,k),:) - vertex(jx,:);
                x2o = vertex(H(jx,1),:) - vertex(jx,:);
                x1n = x1o;
                x2n = x2o;
                if H(jx,k) == ix
                    x1n = x1n + dr;
                end
                if H(jx,1) == ix
                    x2n = x2n + dr;
                end
                if x1o(1,1) > Lx/2.0
                    x1o(1,1)=x1o(1,1)-Lx;
                elseif x1o(1,1) <= -Lx/2.0
                    x1o(1,1)=x1o(1,1)+Lx;
                end
                if x1o(1,2) > Ly/2.0
                    x1o(1,2)=x1o(1,2)-Ly;
                elseif x1o(1,2) <= -Ly/2.0
                    x1o(1,2)=x1o(1,2)+Ly;
                end
                if x1n(1,1) > Lx/2.0
                    x1n(1,1)=x1n(1,1)-Lx;
                elseif x1n(1,1) <= -Lx/2.0
                    x1n(1,1)=x1n(1,1)+Lx;
                end
                if x1n(1,2) > Ly/2.0
                    x1n(1,2)=x1n(1,2)-Ly;
                elseif x1n(1,2) <= -Ly/2.0
                    x1n(1,2)=x1n(1,2)+Ly;
                end
                if x2o(1,1) > Lx/2.0
                    x2o(1,1)=x2o(1,1)-Lx;
                elseif x2o(1,1) <= -Lx/2.0
                    x2o(1,1)=x2o(1,1)+Lx;
                end
                if x2o(1,2) > Ly/2.0
                    x2o(1,2)=x2o(1,2)-Ly;
                elseif x2o(1,2) <= -Ly/2.0
                    x2o(1,2)=x2o(1,2)+Ly;
                end
                if x2n(1,1) > Lx/2.0
                    x2n(1,1)=x2n(1,1)-Lx;
                elseif x2n(1,1) <= -Lx/2.0
                    x2n(1,1)=x2n(1,1)+Lx;
                end
                if x2n(1,2) > Ly/2.0
                    x2n(1,2)=x2n(1,2)-Ly;
                elseif x2n(1,2) <= -Ly/2.0
                    x2n(1,2)=x2n(1,2)+Ly;
                end
                
                CosTheta_o = max(min(dot(x1o,x2o)/(norm(x1o)*norm(x2o)),1),-1);
                CosTheta_n = max(min(dot(x1n,x2n)/(norm(x1n)*norm(x2n)),1),-1);
                theta_o(k,1) = acos(CosTheta_o);
                theta_n(k,1) = acos(CosTheta_n);
            end
        end
      
        if abs(theta_o(1,1) + theta_o(2,1) - theta_o(3,1)) < 1e-8
            theta_o(3,1) = 2*pi - theta_o(3,1);
        end
        if abs(theta_o(1,1) + theta_o(3,1) - theta_o(2,1)) < 1e-8
            theta_o(2,1) = 2*pi - theta_o(2,1);
        end
        if abs(theta_o(2,1) + theta_o(3,1) - theta_o(1,1)) < 1e-8
            theta_o(1,1) = 2*pi - theta_o(1,1);
        end
    
        if abs(theta_n(1,1) + theta_n(2,1) - theta_n(3,1)) < 1e-8
            theta_n(3,1) = 2*pi - theta_n(3,1);
        end
        if abs(theta_n(1,1) + theta_n(3,1) - theta_n(2,1)) < 1e-8
            theta_n(2,1) = 2*pi - theta_n(2,1);
        end
        if abs(theta_n(2,1) + theta_n(3,1) - theta_n(1,1)) < 1e-8
            theta_n(1,1) = 2*pi - theta_n(1,1);
        end
        
        
        for k = 1:ct
            dE_old = dE_old + (theta_o(k,1) - pi*2.0/3.0)^2;
            dE_new = dE_new + (theta_n(k,1) - pi*2.0/3.0)^2;
        end
    end
    
    
    if dE_old > dE_new
        vertex(ix,:) = vn;
        dE(nt,1) = dE_new - dE_old;
        %break
    end
    
    %{
    dE2 = 0;
    dd = 0;
    vertex2 = vertex;
    vertex2(ix,:) = vn;
    for pp = 1:3
        b = vector_p(vertex2(ix,:),vertex2(H(ix,pp),:),Lx,Ly);
        dE2 = (norm(b) - 1) * (norm(b) - 1) + dE2;
    end

    theta_n1 = zeros(3,1);
    for pp = 1:2
        x1n = vector_p(vertex2(ix,:),vertex2(H(ix,pp),:),Lx,Ly);
        x2n = vector_p(vertex2(ix,:),vertex2(H(ix,pp+1),:),Lx,Ly);
        CosTheta_n = max(min(dot(x1n,x2n)/(norm(x1n)*norm(x2n)),1),-1);
        theta_n1(pp,1) = acos(CosTheta_n);
    end
    
    x1n = vector_p(vertex2(ix,:),vertex2(H(ix,3),:),Lx,Ly);
    x2n = vector_p(vertex2(ix,:),vertex2(H(ix,1),:),Lx,Ly);
    CosTheta_n = max(min(dot(x1n,x2n)/(norm(x1n)*norm(x2n)),1),-1);
    theta_n1(3,1) = acos(CosTheta_n);
    
    if abs(theta_n1(1,1) + theta_n1(2,1) - theta_n1(3,1)) < 1e-8
            theta_n1(3,1) = 2*pi - theta_n1(3,1);
    end
    if abs(theta_n1(1,1) + theta_n1(3,1) - theta_n1(2,1)) < 1e-8
            theta_n1(2,1) = 2*pi - theta_n1(2,1);
    end
    if abs(theta_n1(2,1) + theta_n1(3,1) - theta_n1(1,1)) < 1e-8
        theta_n1(1,1) = 2*pi - theta_n1(1,1);
    end
    
    for pp = 1:3
    dE2 = dE2 + (theta_n1(pp,1) - 2/3*pi)^2;
    end


    for pq=1:3
    jx = H(ix,pq);
    theta_n1 = zeros(3,1);
    for pp = 1:2
        x1n = vector_p(vertex2(jx,:),vertex2(H(jx,pp),:),Lx,Ly);
        x2n = vector_p(vertex2(jx,:),vertex2(H(jx,pp+1),:),Lx,Ly);
        CosTheta_n = max(min(dot(x1n,x2n)/(norm(x1n)*norm(x2n)),1),-1);
        theta_n1(pp,1) = acos(CosTheta_n);
    end

    x1n = vector_p(vertex2(jx,:),vertex2(H(jx,3),:),Lx,Ly);
    x2n = vector_p(vertex2(jx,:),vertex2(H(jx,1),:),Lx,Ly);
    CosTheta_n = max(min(dot(x1n,x2n)/(norm(x1n)*norm(x2n)),1),-1);
    theta_n1(3,1) = acos(CosTheta_n);
   
    if abs(theta_n1(1,1) + theta_n1(2,1) - theta_n1(3,1)) < 1e-8
            theta_n1(3,1) = 2*pi - theta_n1(3,1);
    end
    if abs(theta_n1(1,1) + theta_n1(3,1) - theta_n1(2,1)) < 1e-8
            theta_n1(2,1) = 2*pi - theta_n1(2,1);
    end
    if abs(theta_n1(2,1) + theta_n1(3,1) - theta_n1(1,1)) < 1e-8
        theta_n1(1,1) = 2*pi - theta_n1(1,1);
    end
   
    for pp = 1:3
        dE2 = dE2 + (theta_n1(pp,1) - 2/3*pi)^2;
    end
    end
    if abs(dE_new - dE2) > 0.0001
        dd = 1;
        break
    end
    %}
end

A(4:3+N,:) = vertex;
dlmwrite('./vertex.txt', A, 'delimiter','\t','precision',15);

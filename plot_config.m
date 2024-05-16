A=dlmread('./vertex.txt');
%A=dlmread('./vertex_honeycomb.txt');
N=A(1,1);
Lx=A(2,1);
Ly=A(3,2);
vertex=A(4:3+N,:);
H1=dlmread('./connectivity_matrix.txt');
%H1=dlmread('./connectivity_matrix_honeycomb.txt');

for i=1:N
    for j=1:i
        if H1(i,j)==1
            
            dx=abs(vertex(i,1)-vertex(j,1));
            dy=abs(vertex(i,2)-vertex(j,2));
            
            if (dx<Lx/2.0) && (dy<Ly/2.0)
                plot([vertex(i,1) vertex(j,1)],[vertex(i,2) vertex(j,2)],'Color','r','LineWidth',2);hold on;
            elseif (dx<Lx/2.0) && (dy>Ly/2.0)
                if vertex(i,2)<vertex(j,2)
                    plot([vertex(i,1) vertex(j,1)],[vertex(i,2)+Ly vertex(j,2)],'Color','r','LineWidth',2);hold on;
                    plot([vertex(i,1) vertex(j,1)],[vertex(i,2) vertex(j,2)-Ly],'Color','r','LineWidth',2);hold on;
                else
                    plot([vertex(i,1) vertex(j,1)],[vertex(i,2) vertex(j,2)+Ly],'Color','r','LineWidth',2);hold on;
                    plot([vertex(i,1) vertex(j,1)],[vertex(i,2)-Ly vertex(j,2)],'Color','r','LineWidth',2);hold on;
                end
            else
                if (dx>Lx/2.0) && (dy<Ly/2.0)
                    if vertex(i,1)<vertex(j,1)
                        plot([vertex(i,1)+Lx vertex(j,1)],[vertex(i,2) vertex(j,2)],'Color','r','LineWidth',2);hold on;
                        plot([vertex(i,1) vertex(j,1)-Lx],[vertex(i,2) vertex(j,2)],'Color','r','LineWidth',2);hold on;
                    else
                        plot([vertex(i,1)-Lx vertex(j,1)],[vertex(i,2) vertex(j,2)],'Color','r','LineWidth',2);hold on;
                        plot([vertex(i,1) vertex(j,1)+Lx],[vertex(i,2) vertex(j,2)],'Color','r','LineWidth',2);hold on;
                    end
                else
                    if (vertex(i,1)-Lx/2.0) * (vertex(i,2)-Ly/2.0) > 0.0
                        if vertex(i,1)>vertex(j,1)
                            plot([vertex(i,1) vertex(j,1)+Lx],[vertex(i,2) vertex(j,2)+Ly],'Color','r','LineWidth',2);hold on;
                            plot([vertex(i,1) vertex(j,1)+Lx],[vertex(i,2)-Ly vertex(j,2)],'Color','r','LineWidth',2);hold on;
                            plot([vertex(i,1)-Lx vertex(j,1)],[vertex(i,2) vertex(j,2)+Ly],'Color''r','LineWidth',2);hold on;
                            plot([vertex(i,1)-Lx vertex(j,1)],[vertex(i,2)-Ly vertex(j,2)],'Color','r','LineWidth',2);hold on;
                        else
                            plot([vertex(i,1)+Lx vertex(j,1)],[vertex(i,2)+Ly vertex(j,2)],'Color','r','LineWidth',2);hold on;
                            plot([vertex(i,1) vertex(j,1)-Lx],[vertex(i,2)+Ly vertex(j,2)],'Color','r','LineWidth',2);hold on;
                            plot([vertex(i,1)+Lx vertex(j,1)],[vertex(i,2) vertex(j,2)-Ly],'Color','r','LineWidth',2);hold on;
                            plot([vertex(i,1) vertex(j,1)-Lx],[vertex(i,2) vertex(j,2)-Ly],'Color','r','LineWidth',2);hold on;
                        end
                    else
                        if vertex(i,1)>vertex(j,1)
                            plot([vertex(i,1) vertex(j,1)+Lx],[vertex(i,2) vertex(j,2)-Ly],'Color','r','LineWidth',2);hold on;
                            plot([vertex(i,1)-Lx vertex(j,1)],[vertex(i,2) vertex(j,2)-Ly],'Color','r','LineWidth',2);hold on;
                            plot([vertex(i,1)-Lx vertex(j,1)],[vertex(i,2)+Ly vertex(j,2)],'Color','r','LineWidth',2);hold on;
                            plot([vertex(i,1) vertex(j,1)+Lx],[vertex(i,2)+Ly vertex(j,2)],'Color','r','LineWidth',2);hold on;
                        else
                            plot([vertex(i,1)+Lx vertex(j,1)],[vertex(i,2)-Ly vertex(j,2)],'Color','r','LineWidth',2);hold on;
                            plot([vertex(i,1) vertex(j,1)-Lx],[vertex(i,2) vertex(j,2)+Ly],'Color','r','LineWidth',2);hold on;
                            plot([vertex(i,1) vertex(j,1)-Lx],[vertex(i,2)-Ly vertex(j,2)],'Color','r','LineWidth',2);hold on;
                            plot([vertex(i,1)+Lx vertex(j,1)],[vertex(i,2) vertex(j,2)+Ly],'Color','r','LineWidth',2);hold on;
                        end
                    end
                end
            end
        end
    end
end


axis equal;
axis([0 Lx 0 Ly]);
set(gca,'xtick',[],'ytick',[]);
set(gca,'FontSize',15);
plot([0.0 0.0], [0.0 Ly], 'Color', 'k', 'LineWidth',2);
plot([Lx Lx], [0.0 Ly], 'Color', 'k', 'LineWidth',2);
plot([0.0 Lx], [0.0 0.0], 'Color', 'k', 'LineWidth',2);
plot([0.0 Lx], [Ly Ly], 'Color', 'k', 'LineWidth',2);

% E is a collection of tetrahedra filling the 3D mesh.
% That is, a list of quadruples of points (numbered top-to-bottom, 
%   then left-to-right, then front-to-back) which form tetrahedra.
% As appropriate, we record the tetrahedra, the faces, and the edges.
function [edges,faces,tetrahedra]=TetrahedralMesh(l1,l2,l3,dim)
n1=l1+1;
n2=l2+1;
n3=l3+1;

number_of_tetrahedra = l1*l2*l3*6;
tetrahedra=zeros(number_of_tetrahedra,4);

number_of_faces = l1*l2*l3*12 + l1*l2*2 + l1*l3*2 + l2*l3*2;
faces=zeros(number_of_faces,3);

number_of_edges = l1*l2*l3*7 + l1*l2*3 + l1*l3*3 + l2*l3*3 + l1 + l2 + l3;
edges=zeros(number_of_edges,2);

counter1 = 1; % counter for the tetrahedra
counter2 = 1; % counter for the triangular faces
counter3 = 1; % counter for the edges
for i=1:n3-1
    for j=1:n2-1
        for k=1:n1-1
            v = (i-1)*n1*n2+(j-1)*n1+k;
            
            if (dim==2)
                tetrahedra(counter1,:)=[v v+1 v+n1+1 v+n1*n2+n1+1]; % 6 terahedra that start with the same points v
                tetrahedra(counter1+1,:)=[v v+1 v+n1*n2+1 v+n1*n2+n1+1];
                tetrahedra(counter1+2,:)=[v v+n1 v+n1+1 v+n1*n2+n1+1];
                tetrahedra(counter1+3,:)=[v v+n1 v+n1*n2+n1 v+n1*n2+n1+1];
                tetrahedra(counter1+4,:)=[v v+n1*n2 v+n1*n2+1 v+n1*n2+n1+1];
                tetrahedra(counter1+5,:)=[v v+n1*n2 v+n1*n2+n1 v+n1*n2+n1+1];
                counter1 = counter1 + 6;
            end
            
            if (dim==1) || (dim==2)
                faces(counter2,:)=[v v+1 v+n1+1]; % 12 faces that start with the same point v
                faces(counter2+1,:)=[v v+1 v+n1*n2+1];
                faces(counter2+2,:)=[v v+1 v+n1*n2+n1+1];
                faces(counter2+3,:)=[v v+n1 v+n1+1];
                faces(counter2+4,:)=[v v+n1 v+n1*n2+n1];
                faces(counter2+5,:)=[v v+n1 v+n1*n2+n1+1];
                faces(counter2+6,:)=[v v+n1+1 v+n1*n2+n1+1];
                faces(counter2+7,:)=[v v+n1*n2 v+n1*n2+1];
                faces(counter2+8,:)=[v v+n1*n2 v+n1*n2+n1];
                faces(counter2+9,:)=[v v+n1*n2 v+n1*n2+n1+1];
                faces(counter2+10,:)=[v v+n1*n2+1 v+n1*n2+n1+1];
                faces(counter2+11,:)=[v v+n1*n2+n1 v+n1*n2+n1+1];
                counter2 = counter2 + 12;
            end
            
            edges(counter3,:)=[v v+1]; % 7 edges that start with the same point v
            edges(counter3+1,:)=[v v+n1];
            edges(counter3+2,:)=[v v+n1+1];
            edges(counter3+3,:)=[v v+n1*n2];
            edges(counter3+4,:)=[v v+n1*n2+1];
            edges(counter3+5,:)=[v v+n1*n2+n1];
            edges(counter3+6,:)=[v v+n1*n2+n1+1];
            counter3 = counter3 + 7;
        end
        % In the last page
        k = n1;
        v = (i-1)*n1*n2+(j-1)*n1+k;
        
        if (dim==1) || (dim==2)
            faces(counter2,:)=[v v+n1 v+n1*n2+n1];
            faces(counter2+1,:)=[v v+n1*n2 v+n1*n2+n1];
            counter2 = counter2 + 2;
        end
        
        edges(counter3,:)=[v v+n1];
        edges(counter3+1,:)=[v v+n1*n2];
        edges(counter3+2,:)=[v v+n1*n2+n1];
        counter3 = counter3 + 3;
    end
    % In the last column
    j = n2;
    for k=1:n1-1
        v = (i-1)*n1*n2+(j-1)*n1+k;
        
        if (dim==1) || (dim==2)
            faces(counter2,:)=[v v+1 v+n1*n2+1];
            faces(counter2+1,:)=[v v+n1*n2 v+n1*n2+1];
            counter2 = counter2 + 2;
        end
        
        edges(counter3,:)=[v v+1];
        edges(counter3+1,:)=[v v+n1*n2];
        edges(counter3+2,:)=[v v+n1*n2+1];
        counter3 = counter3 + 3;
    end
    
    k = n1;
    v = (i-1)*n1*n2+(j-1)*n1+k;
    
    edges(counter3,:)=[v v+n1*n2];
    counter3 = counter3 + 1;
end
% In the last row
i = n3;
for j=1:n2-1
    for k=1:n1-1
        v = (i-1)*n1*n2+(j-1)*n1+k;
        
        if (dim==1) || (dim==2)
            faces(counter2,:)=[v v+1 v+n1+1];
            faces(counter2+1,:)=[v v+n1 v+n1+1];
            counter2 = counter2 + 2;
        end
        
        edges(counter3,:)=[v v+1];
        edges(counter3+1,:)=[v v+n1];
        edges(counter3+2,:)=[v v+n1+1];
        counter3 = counter3 + 3;
    end
    k = n1;
    v = (i-1)*n1*n2+(j-1)*n1+k;
    
    edges(counter3,:)=[v v+n1];
    counter3 = counter3 + 1;
end
j = n2;
for k=1:n1-1
    v = (i-1)*n1*n2+(j-1)*n1+k;
    
    edges(counter3,:)=[v v+1];
    counter3 = counter3 + 1;
end

end
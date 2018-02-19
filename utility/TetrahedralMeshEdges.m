% E is a collection of tetrahedra filling the 3D mesh.
% That is, a list of quadruples of points (numbered top-to-bottom, 
%   then left-to-right, then front-to-back) which form tetrahedra.
% As appropriate, we record the tetrahedra, the faces, and the edges.
function edges=TetrahedralMeshEdges(n1,n2,n3,dim)
l1=n1-1;
l2=n2-1;
l3=n3-1;


number_of_edges = l1*l2*l3*7 + l1*l2*3 + l1*l3*3 + l2*l3*3 + l1 + l2 + l3;
edges=zeros(number_of_edges,2);

counter = 1; % counter for the edges
for i=1:n3-1
    for j=1:n2-1
        for k=1:n1-1
            v = (i-1)*n1*n2+(j-1)*n1+k;
            
            edges(counter,:)=[v v+1]; % 7 edges that start with the same point v
            edges(counter+1,:)=[v v+n1];
            edges(counter+2,:)=[v v+n1+1];
            edges(counter+3,:)=[v v+n1*n2];
            edges(counter+4,:)=[v v+n1*n2+1];
            edges(counter+5,:)=[v v+n1*n2+n1];
            edges(counter+6,:)=[v v+n1*n2+n1+1];
            counter = counter + 7;
        end
        % In the last page
        k = n1;
        v = (i-1)*n1*n2+(j-1)*n1+k;
        
        edges(counter,:)=[v v+n1];
        edges(counter+1,:)=[v v+n1*n2];
        edges(counter+2,:)=[v v+n1*n2+n1];
        counter = counter + 3;
    end
    % In the last column
    j = n2;
    for k=1:n1-1
        v = (i-1)*n1*n2+(j-1)*n1+k;

        edges(counter,:)=[v v+1];
        edges(counter+1,:)=[v v+n1*n2];
        edges(counter+2,:)=[v v+n1*n2+1];
        counter = counter + 3;
    end
    % In the last page and the last column
    k = n1;
    v = (i-1)*n1*n2+(j-1)*n1+k;
    
    edges(counter,:)=[v v+n1*n2];
    counter = counter + 1;
end
% In the last row
i = n3;
for j=1:n2-1
    for k=1:n1-1
        v = (i-1)*n1*n2+(j-1)*n1+k;

        edges(counter,:)=[v v+1];
        edges(counter+1,:)=[v v+n1];
        edges(counter+2,:)=[v v+n1+1];
        counter = counter + 3;
    end
    % In the last column and the last row
    k = n1;
    v = (i-1)*n1*n2+(j-1)*n1+k;
    
    edges(counter,:)=[v v+n1];
    counter = counter + 1;
end
% In the last row and last page
j = n2;
for k=1:n1-1
    v = (i-1)*n1*n2+(j-1)*n1+k;
    
    edges(counter,:)=[v v+1];
    counter = counter + 1;
end

end
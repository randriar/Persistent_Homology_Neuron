% method that computes the Maximum Euclidean Distance of a neuron
% Parameter: cell- the name of the cell

% We get the Euclidean Distance by computing the straight line distance
%   from the soma (root) to the node 
%   => sqrt((x-somax)^2+(y-somay)^2+(z-somaz)^2)
function f=NeuronRadius(cell) % Maximum euclidean distance from the soma
soma=cell(1,:); % soma (the first point of the cell)
x_coords=cell(:,1); % list of x-coordinates of the cell
y_coords=cell(:,2); % list of y-coordinates of the cell
z_coords=cell(:,3); % list of z-coordinates of the cell
x=x_coords-soma(:,1); 
y=y_coords-soma(:,2);
z=z_coords-soma(:,3);
distances=sqrt(x.^2+y.^2+z.^2); % Compute the Euclidean distances between each node and the soma
max_distance=max(distances); % take the maximum distance
f=max_distance;

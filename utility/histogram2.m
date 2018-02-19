% ######## ************ #######
% PARAMETERS:
%   cell: x,y and z corrdinates of the cell
%   rad: list of radii
%   gridx, gridy, gridz

% OUTPUT:
%   The kernel density of the cell presented in a 3D-matrix 
% ######## ************ #######

function f=histogram2(cell,rad,gridx,gridy,gridz)

cellx= cell(:,1); % the x-coordinates of the cell
celly= cell(:,2); % the y-coordinates of the cell
cellz= cell(:,3); % the z-coordinates of the cell

nb=1; % initialize a variable that keeps track of the incrementation when adding new points

% Create a new array and store the new points created from the radius of
%   each points in it.
%   28 points are created from the radius for each point in the cell
new_pts= zeros(26*size(rad,1),3); 
for l=1:1:(size(rad,1))
    
    new_pts(nb  ,:)= [ cellx(l)+rad(l)            celly(l)                   cellz(l)                  ];
    new_pts(nb+1,:)= [ cellx(l)-rad(l)            celly(l)                   cellz(l)                  ];
    new_pts(nb+2,:)= [ cellx(l)                   celly(l)+rad(l)            cellz(l)                  ];
    new_pts(nb+3,:)= [ cellx(l)                   celly(l)-rad(l)            cellz(l)                  ];
    new_pts(nb+4,:)= [ cellx(l)                   celly(l)                   cellz(l)+rad(l)           ];
    new_pts(nb+5,:)= [ cellx(l)                   celly(l)                   cellz(l)-rad(l)           ];
    
    
    new_pts(nb+6,:)= [ cellx(l)+(rad(l)/sqrt(2))  celly(l)+(rad(l)/sqrt(2))  cellz(l)                  ];
    new_pts(nb+7,:)= [ cellx(l)+(rad(l)/sqrt(2))  celly(l)-(rad(l)/sqrt(2))  cellz(l)                  ];
    new_pts(nb+8,:)= [ cellx(l)-(rad(l)/sqrt(2))  celly(l)+(rad(l)/sqrt(2))  cellz(l)                  ];
    new_pts(nb+9,:)= [ cellx(l)-(rad(l)/sqrt(2))  celly(l)-(rad(l)/sqrt(2))  cellz(l)                  ];
    
    new_pts(nb+10,:)=[ cellx(l)+(rad(l)/sqrt(2))  celly(l)                   cellz(l)+(rad(l)/sqrt(2)) ];
    new_pts(nb+11,:)=[ cellx(l)+(rad(l)/sqrt(2))  celly(l)                   cellz(l)-(rad(l)/sqrt(2)) ];    
    new_pts(nb+12,:)=[ cellx(l)-(rad(l)/sqrt(2))  celly(l)                   cellz(l)+(rad(l)/sqrt(2)) ];
    new_pts(nb+13,:)=[ cellx(l)-(rad(l)/sqrt(2))  celly(l)                   cellz(l)-(rad(l)/sqrt(2)) ];
    
    new_pts(nb+14,:)=[ cellx(l)                   celly(l)+(rad(l)/sqrt(2))  cellz(l)+(rad(l)/sqrt(2)) ];
    new_pts(nb+15,:)=[ cellx(l)                   celly(l)+(rad(l)/sqrt(2))  cellz(l)-(rad(l)/sqrt(2)) ];
    new_pts(nb+16,:)=[ cellx(l)                   celly(l)-(rad(l)/sqrt(2))  cellz(l)+(rad(l)/sqrt(2)) ];
    new_pts(nb+17,:)=[ cellx(l)                   celly(l)-(rad(l)/sqrt(2))  cellz(l)-(rad(l)/sqrt(2)) ];
    
    new_pts(nb+18,:)=[ cellx(l)+(rad(l)/sqrt(3))  celly(l)+(rad(l)/sqrt(3))  cellz(l)+(rad(l)/sqrt(3)) ];
    new_pts(nb+19,:)=[ cellx(l)+(rad(l)/sqrt(3))  celly(l)+(rad(l)/sqrt(3))  cellz(l)-(rad(l)/sqrt(3)) ];
    new_pts(nb+20,:)=[ cellx(l)+(rad(l)/sqrt(3))  celly(l)-(rad(l)/sqrt(3))  cellz(l)+(rad(l)/sqrt(3)) ];
    new_pts(nb+21,:)=[ cellx(l)+(rad(l)/sqrt(3))  celly(l)-(rad(l)/sqrt(3))  cellz(l)-(rad(l)/sqrt(3)) ];
    new_pts(nb+22,:)=[ cellx(l)-(rad(l)/sqrt(3))  celly(l)+(rad(l)/sqrt(3))  cellz(l)+(rad(l)/sqrt(3)) ];
    new_pts(nb+23,:)=[ cellx(l)-(rad(l)/sqrt(3))  celly(l)+(rad(l)/sqrt(3))  cellz(l)-(rad(l)/sqrt(3)) ];
    new_pts(nb+24,:)=[ cellx(l)-(rad(l)/sqrt(3))  celly(l)-(rad(l)/sqrt(3))  cellz(l)+(rad(l)/sqrt(3)) ];
    new_pts(nb+25,:)=[ cellx(l)-(rad(l)/sqrt(3))  celly(l)-(rad(l)/sqrt(3))  cellz(l)-(rad(l)/sqrt(3)) ];
    
    nb=nb+26; % increment the index 
end

cell= cat(1,cell,new_pts); % concatenate the two matrices (old cell & new_pts)

cellx= cell(:,1); % the new x-coordinates of the cell
celly= cell(:,2); % the new y-coordinates of the cell
cellz= cell(:,3); % the new z-coordinates of the cell

% Some basic size values and parameters.
[n,~] = size(cell); % size of the sample
nx = length(gridx); % length of the grid x
ny = length(gridy); % length of the grid y
nz = length(gridz); % length of the grid z

max_boundx = gridx(nx); % maximum bound in the x-coordinate
max_boundy = gridy(ny); % max bound in the y-coordinate
max_boundz = gridz(nz); % max bound in the z-coordinate
min_boundx = gridx(1);  % minimum bound in the x-coordinate
min_boundy = gridy(1);  % min bound in the y-coordinate
min_boundz = gridz(1);  % min bound in the z-coordinate

stepx= nx-1; % number of steps along the x-axis
stepy= ny-1; % number of steps along the y-axis
stepz= nz-1; % number of steps along the z-axis
bin_widthx= (max_boundx-min_boundx)/stepx; % length of each bin width along the x-axis
bin_widthy= (max_boundy-min_boundy)/stepy; % length of each bin width along the y-axis
bin_widthz= (max_boundz-min_boundz)/stepz; % length of each bin width along the z-axis

density= zeros(nx,ny,nz); % create a matrix of zeroes to store the density

% Computing the density estimator for each bin in the plan
%   and store them in the appropriate index
% initializing the pivot for z-axis
for k=1:1:stepz 
    z=gridz(k);
    start_binz=z-(bin_widthz/2);
    stop_binz =z+(bin_widthz/2);
    % initializing the pivot for y-axis
    for j=1:1:stepy
        y=gridy(j);
        start_biny=y-(bin_widthy/2);
        stop_biny =y+(bin_widthy/2);
        % initializing the pivot for x-axis
        for i=1:1:stepx
            x= gridx(i);
            start_binx= x - (bin_widthx /2);
            stop_binx = x + (bin_widthx /2);
            indexes= find(cellx>start_binx & cellx<stop_binx & celly>start_biny & celly<stop_biny & cellz>start_binz & cellz<stop_binz);
            nbOfPoints= length(indexes); % number of points in the bin width
            if nbOfPoints >= 1
                density(i,j,k)= nbOfPoints;
            end
        end
    end  
end

f = density; % Return the 3D matrix that represents the density of the cell

% ######## ************ #######
% PARAMETERS:
%   cell
%   r 
%   gridx, gridy, gridz

% OUTPUT:
%   The kernel density of the cell presented in a 3D-matrix 
% ######## ************ #######

function f=histogram(cell,r,gridx,gridy,gridz)

% Some basic size values and parameters.
[n,~] = size(cell); % size of the sample
nx = length(gridx); % length of the grid x
ny = length(gridy); % length of the grid y
nz = length(gridz); % length of the grid z

max_boundx=gridx(nx); % maximum bound in the x-coordinate
max_boundy=gridy(ny); % max bound in the y-coordinate
max_boundz=gridz(nz); % max bound in the z-coordinate
min_boundx=gridx(1);             % minimum bound in the x-coordinate
min_boundy=gridy(1);             % min bound in the y-coordinate
min_boundz=gridz(1);             % min bound in the z-coordinate

stepx=nx-1; % number of steps along the x-axis
stepy=ny-1; % number of steps along the y-axis
stepz=nz-1; % number of steps along the z-axis

cellx=cell(:,1); % the x-coordinates of the cell
celly=cell(:,2); % the y-coordinates of the cell
cellz=cell(:,3); % the z-coordinates of the cell

bin_widthx=(max_boundx-min_boundx)/stepx; % length of each bin width along the x-axis
bin_widthy=(max_boundy-min_boundy)/stepy; % length of each bin width along the y-axis
bin_widthz=(max_boundz-min_boundz)/stepz; % length of each bin width along the z-axis

volume=1; %bin_widthx*bin_widthy*bin_widthz;

density = zeros(nx,ny,nz); % create a matrix of zeroes
% Computing the density estimator for each bin in the plan
%   and store them in the appropriate index
% initializing the pivot for z-axis
for k=1:1:stepz 
    z=gridz(k);
    start_binz=z-(bin_widthz/2);
    stop_binz=z+(bin_widthz/2);
    % initializing the pivot for y-axis
    for j=1:1:stepy
        y=gridy(j);
        start_biny=y-(bin_widthy/2);
        stop_biny=y+(bin_widthy/2);
        % initializing the pivot for x-axis
        for i=1:1:stepx
            x=gridx(i);
            start_binx=x-(bin_widthx/2);
            stop_binx=x+(bin_widthx/2);
            indexes=find(cellx>start_binx & cellx<stop_binx & celly>start_biny & celly<stop_biny & cellz>start_binz & cellz<stop_binz);
            nbOfPoints=length(indexes); % number of points in the bin width
            if nbOfPoints >= 1
                density(i,j,k)=sum(r(indexes)); % Computing the density by using the histogram method
            end
        end
    end  
end

f = density; % Return the 3D matrix that represents the density of the cell

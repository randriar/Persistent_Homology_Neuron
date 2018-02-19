function f=sphere_histogram(rad,gridx,gridy,gridz)
nx = length(gridx); % length of the grid x
ny = length(gridy); % length of the grid y
nz = length(gridz); % length of the grid z

max_boundx=gridx(nx); % maximum bound in the x-coordinate
max_boundy=gridy(ny); % max bound in the y-coordinate
max_boundz=gridz(nz); % max bound in the z-coordinate
min_boundx=gridx(1);  % minimum bound in the x-coordinate
min_boundy=gridy(1);  % min bound in the y-coordinate
min_boundz=gridz(1);  % min bound in the z-coordinate

stepx=nx-1; % number of steps along the x-axis
stepy=ny-1; % number of steps along the y-axis
stepz=nz-1; % number of steps along the z-axis

bin_widthx=(max_boundx-min_boundx)/stepx; % length of each bin width along the x-axis
bin_widthy=(max_boundy-min_boundy)/stepy; % length of each bin width along the y-axis
bin_widthz=(max_boundz-min_boundz)/stepz; % length of each bin width along the z-axis

density = zeros(nx,ny,nz); % create a matrix of zeroes

for k=1:1:nz
    c=gridz(k);
    start_binz=c-(bin_widthz/2);
    stop_binz=c+(bin_widthz/2);
    % initializing the pivot for y-axis
    for j=1:1:ny
        b=gridy(j);
        start_biny=b-(bin_widthy/2);
        stop_biny=b+(bin_widthy/2);
        % initializing the pivot for x-axis
        for i=1:1:nx
            a=gridx(i);
            start_binx=a-(bin_widthx/2);
            stop_binx=a+(bin_widthx/2);
            
            min_x = min(abs(start_binx),abs(stop_binx));
            max_x = max(abs(start_binx),abs(stop_binx));
            min_y = min(abs(start_biny),abs(stop_biny));
            max_y = max(abs(start_biny),abs(stop_biny));
            min_z = min(abs(start_binz),abs(stop_binz));
            max_z = max(abs(start_binz),abs(stop_binz));
            
            min_bin=min_x^2+min_y^2+min_z^2;
            max_bin=max_x^2+max_y^2+max_z^2;
            if( (rad^2)>=min_bin && (rad^2)<max_bin )
                density(i,j,k)=1;
            end
        end
    end
end

f=density;
% source code: Process any 3D cell and compute its Persistent Homology
%   methods used: LoadNeuron
%                 NeuronRadius - Euclidian Distance
%                 TetrahedralMesh
%                 sphere_histogram - the KDE of the sphere
%                 histogram2 - KDE of the the input cell
%                 ShollAnalysis
%                 TraditionalMorphology
%
% author: Prof. Carl Hammarsten & Rindra Randriamanantena
% version: February 12, 2018
% 

% This script prepares the javaplex library for use.
clc; clear all; close all;
javaaddpath('./utility/javaplex.jar');
import edu.stanford.math.plex4.*;
javaaddpath('./utility/plex-viewer.jar');
import edu.stanford.math.plex_viewer.*;
cd './utility';
addpath(pwd);
cd '..';

% Set some parameters.
resolutionx = 5;
resolutiony = 5;
resolutionz = 5;
radii_list = [25;50;75;100;125;150;200;250;300;400;450;475;500;525;550];
dim = 0;

disp('START...');

% Create directories where the data will be saved
if exist('OUTPUTS','file')==7   % If the 'OUTPUTS' folder already exists, delete it
    rmdir('OUTPUTS','s');
end
mkdir('OUTPUTS');               % create a directory where the data about each cells will be saved
mkdir('OUTPUTS','Parameters');  % create a directory to save the parameters
% mkdir('OUTPUTS','SphericalShells'); % create a directory to save the spherical shells
mkdir('OUTPUTS','CELLS');           % Create a folder called 'CELLS' where the data about each cell will be stored
save(strcat('OUTPUTS/Parameters/resolutionx'),'resolutionx');
save(strcat('OUTPUTS/Parameters/resolutiony'),'resolutiony');
save(strcat('OUTPUTS/Parameters/resolutionz'),'resolutionz');
save(strcat('OUTPUTS/Parameters/radii'),'radii_list');
save(strcat('OUTPUTS/Parameters/dimension'),'dim');
% mkdir('OUTPUTS/SphericalShells','SphericalShell_List');    % Create a directory to save the lists of spherical shells
% mkdir('OUTPUTS/SphericalShells','SphericalShell_Matrix');  % Create a directory to save the matrices of spherical shells

disp('Processing the Cells');
file=dir('INPUTS/*.swc'); % Take the files with '.swc' extension
for k=1:length(file) % For every cell
    disp(strcat(' |--Cell #',num2str(k)));
    
    filename=file(k).name;
    cell0=LoadNeuron(strcat('INPUTS/',filename)); % load the neuron in the 'INPUTS' folder
    % cell1=cell0;
    % Making directories where to save the data
    filename=strsplit(filename,'.swc'); % name of the cell without extension
    filename=filename{1}; % filename without extension
    mkdir('OUTPUTS/CELLS',strcat('Cell#',num2str(k),'_',filename)); % Create a directory where the Cell's information(pictures,etc...) will be stored
    mkdir(strcat('OUTPUTS/CELLS/Cell#',num2str(k),'_',filename),'Sphere*KDE_Matrix');   % Create a directory to save the product of the sphere and KDE (matrices)
    mkdir(strcat('OUTPUTS/CELLS/Cell#',num2str(k),'_',filename),'Sphere*KDE_List');     % Create a directory to save the product of the spherical shells and KDE (lists)
    mkdir(strcat('OUTPUTS/CELLS/Cell#',num2str(k),'_',filename),'Homology');            % Create a directory to save the Homology(barcodes) of the cell     
    mkdir(strcat('OUTPUTS/CELLS/Cell#',num2str(k),'_',filename,'/Homology'),'Barcodes');% Create a directory to save the barcodes
    mkdir(strcat('OUTPUTS/CELLS/Cell#',num2str(k),'_',filename),'Images');              % Create a directory to save images
    
    id_list=cell0(:,1);
    parent_list=cell0(:,7); % setting the parents to the 4th column of the cell
    r=cell0(:,6);
    cell0(:,7)=[]; % emptying the 7th column
    cell0(:,6)=[]; % emptying the 6th column
    cell0(:,2)=[]; % emptying the 2th column
    cell0(:,1)=[]; % emptying the 4th column
    save(strcat('OUTPUTS/CELLS/','Cell#',num2str(k),'_',filename,'/XYZ_Coordinates_Raw'),'cell0'); % Saving the cell
    
    % Compute the maximum Euclidean distance
    % For more information, see 'NeuronRadius.m'
    max_Euclidean_distance=NeuronRadius(cell0); 

    % Recenter at the "soma"
    disp(' |   |--Recentering Cell');
    
    % Shift so that the "physical center of the cell" - AKA the "soma" is
    %   at (0,0,0). Our standard Cell format should include the soma
    %   coordinates as the first coordinate-triple.
    cell = cell0-repmat(cell0(1,:),length(cell0),1); % recenter the cell 
    
    % Save a scatter plot of the cells points comprising this neuron.
    fig0 = figure('visible', 'off'); % turn off the visibility of the figure
    scatter3(cell(:,1),cell(:,2),cell(:,3),'.'); % plot the original cell in a 3D plan 
    print(strcat('OUTPUTS/CELLS/','Cell#',num2str(k),'_',filename,'/Images/Cell_Plot'),'-dpdf'); % save the 3D image of the cell
    close(fig0);
    
    % Create viewing window. Goes one extra "resolution sized" step in
    % each direction in every dimension.
    clear U;
    lowerboundx = (floor(min(cell(:,1))/resolutionx)-1)*resolutionx;
    upperboundx = (ceil(max(cell(:,1))/resolutionx)+1)*resolutionx;
    lowerboundy = (floor(min(cell(:,2))/resolutiony)-1)*resolutiony;
    upperboundy = (ceil(max(cell(:,2))/resolutiony)+1)*resolutiony;
    lowerboundz = (floor(min(cell(:,3))/resolutionz)-1)*resolutionz;
    upperboundz = (ceil(max(cell(:,3))/resolutionz)+1)*resolutionz;
    
    sigma=sqrt((resolutionx/2)^2+(resolutiony/2)^2+(resolutionz/2)^2); % Standard deviation used to determine "thickness" of spherical shells.
    
    hgrid= lowerboundx:resolutionx:upperboundx;          % Discrete points along x axis.
    vgrid= lowerboundy:resolutiony:upperboundy;          % Discrete points along y axis.
    dgrid= lowerboundz:resolutionz:upperboundz;          % Discrete points along z axis.
    
    [xsample,ysample,zsample]=meshgrid(hgrid,vgrid,dgrid);       % Mesh in 3D space (copies of hgrid, vgrid and dgrid).
    
    U(:,1)= reshape(xsample,size(xsample,1)*size(xsample,2)*size(xsample,3),1);  % Re-format, the mesh as a list of the
    U(:,2)= reshape(ysample,size(ysample,1)*size(ysample,2)*size(ysample,3),1);  %  coordinate-triples 
    U(:,3)= reshape(zsample,size(zsample,1)*size(zsample,2)*size(zsample,3),1);  %  for the vertices 
    
    n1= length(hgrid);n2=length(vgrid);n3=length(dgrid); % The number of points in each directional axis.
    
    disp(' |   |--Creating Adjacency Mesh');
    edges= TetrahedralMeshEdges(n1,n2,n3,dim);
    number_of_edges= length(edges);
    
    % Compute the Gaussian density estimator. See 'KernelDensityEstimator.m' for details.
    disp(' |   |--Computing KDE');
    KDE=histogram2(cell,r,hgrid,vgrid,dgrid); % computing kde by using the histogram method

    save(strcat('OUTPUTS/CELLS/Cell#',num2str(k),'_',filename,'/Cell_KDE'),'KDE'); % saving the KDE
    
    disp(' |   |--Computing Spherical KDEs and Barcodes');
    for t=1:length(radii_list) % for each spherical shells
        disp(strcat(' |   |   |--Sphere #',num2str(t)))
        rad = radii_list(t);

        % If current spherical radius is less than the neuron radius, we check for
        %   homology.
        if(rad<=max_Euclidean_distance)
            fig0 = figure;
            [x,y,z] = sphere(17);
            x = rad*x;
            y = rad*y;
            z = rad*z;
            hold on;
            scatter3(cell(:,1),cell(:,2),cell(:,3),'.');
            surf(x,y,z);
            alpha .05;
            print(strcat('OUTPUTS/CELLS/','Cell#',num2str(k),'_',filename,'/Images/2D_Cell+Sphere',num2str(rad)),'-dpdf');
            close(fig0);
            
            disp(' |   |   |   |--Computing Spherical KDE');
            
            % Compute the spherical shell kernel.
            SphericalShell=sphere_histogram(rad,hgrid,vgrid,dgrid); 

            % Multiply the cell KDE by the spherical shell kernel.
            product=reshape(reshape(KDE,size(SphericalShell)).*SphericalShell,size(KDE)); 
            ZZ=-reshape(product,n1*n2*n3,1);                                                                                                     % Negative function - switch superlevel set to sublevel set
            save(strcat('OUTPUTS/CELLS/Cell#',num2str(k),'_',filename,'/Sphere*KDE_Matrix/SphericalShell*KDE_Matrix_', num2str(t)),'product');   % save the product as a matrix
            save(strcat('OUTPUTS/CELLS/Cell#',num2str(k),'_',filename,'/Sphere*KDE_List/SphericalShell*KDE_List_', num2str(t)),'ZZ');            % save the product as a list

            % Compute persistent barcodes for each spherical shell.
            %   Initialize an explicit chain complex.
            disp(' |   |   |   |--Creating Simplex');
            stream = api.Plex4.createExplicitSimplexStream(2); 
            
            % Add vertices to our explicit chain complex. 
            %   Each vertex is given a 'birth time' according to the
            %   KDE*sphere value at this point.
            for i=1:length(U) % adding vertices
                if ZZ(i) < 0
                    stream.addVertex(i,ZZ(i));
                end
            end
            
            % Add edges given 2 vertices which are adjacent in our 3D mesh.
            %   Each edge is given a 'birth time' according to the larger
            %   KDE*sphere at any end vertex.
            for i=1:number_of_edges
                if (ZZ(edges(i,1)) < 0) || (ZZ(edges(i,2)) < 0)
                    stream.addElement([edges(i,1) edges(i,2)],max([ZZ(edges(i,1)),ZZ(edges(i,2))]));
                end
            end
    
            % Finalize the chain complex.
            stream.finalizeStream();
        
            disp(' |   |   |   |--Computing Barcode');
            % Compute homology in dimension "dim", and with coefficients in Z_2.
            persistence = api.Plex4.getModularSimplicialAlgorithm((dim+1), 2);
            % Compute persistence intervals.
            intervals = persistence.computeIntervals(stream);
            % Some options for displaying the barcodes.
            options.min_dimension=0;
            options.max_dimension=dim;
            options.min_filtration_value=min(ZZ);
            options.max_filtration_value=0;
            plot_barcodes(intervals,options); % Actually display the barcodes.
            
            print(strcat('OUTPUTS\CELLS\','Cell#',num2str(k),'_',filename,'\Homology\Barcodes\','Barcode_', num2str(rad)),'-dpdf');   % saving the barcodes in the corresponding folder of the cooresponding cell 
            close;
            Interval0{k}{t}=intervals.getIntervalsAtDimension(0); % Record the 0-dim intervals for pairwise comparison.
            
        % If current radius is larger than neuron radius, we set all homology to null.    
        %   i.e. There is no cell to look at in this particular region.
        else
            disp(' |   |   |   |--Radius Larger Than Cell');
            Interval0{k}{t}=java.util.ArrayList;
        end
    end
    
    % Save the homology structures.
    save(strcat('OUTPUTS/CELLS/','Cell#',num2str(k),'_',filename,'/Homology/Interval0'),'Interval0') % Save the 0-dimensional homology in Interval0.
    
    % Compute traditional Morphology measurements.
    ReturnPath = cd(strcat('./OUTPUTS/CELLS/Cell#',num2str(k),'_',filename));
    disp(' |   |--Computing Fuzzy Sholl Analysis');
    ShollAnalysis(Interval0{k},radii_list);
    % save(strcat('OUTPUTS/CELLS/','Cell#',num2str(k),'_',filename,'/Homology/Interval0'),'Interval0') % Save the Sholl Analysis
    %print(strcat('OUTPUTS\CELLS\','Cell#',num2str(k),'_',filename,'\ShollAnalysis\','Sholl'),'-dpdf');
    disp(' |   |--Computing Traditional Morphology Indicators');
    TraditionalMorphology(id_list,cell0,parent_list);
    cd(ReturnPath);
end

disp('FINISHED!');
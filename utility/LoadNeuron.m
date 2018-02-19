%                   ###### ********* ######

%       This method load cells in swc format so that we can use them in the main part of the program 
%       Parameter: cell- the name of the cell

%                   ###### ********* ######

function f=LoadNeuron(cell)
% Change the extension to '.txt'
filename=strsplit(cell,'.swc'); % split the filename
filenameWithoutExtension=filename{1}; % filename without the extension
newFile=strcat(filenameWithoutExtension,'.txt'); % new filename with the new extension 'txt'
copyfile(cell,newFile); % make a copy of the file and change the extension to txt 
fid=fopen(newFile); % open the new file
C= textscan(fid, '%s','delimiter', '\n'); % scan the new file and remove the 3 first lines
delete(newFile); % delete the copy 
D = C{1}; % the part the need
stringCell=string(D); % Convert Cell to String
index=-1; % index of the line where the cell starts
i=1;
while index==-1 % find index
  L=strsplit(stringCell(i),{' ',','}); % split each line by using ' ' and ',' as delimiters
  if(double(L(1))==1 || double(L(1))==0) % find the index where the cell starts
      index=i;
      break; % stop the loop when the index is found
  end
  i=i+1;
end 
cell = zeros(length(stringCell)-index+1,7); % create a cell of 7 columns with 0s elements
for j=1:length(stringCell)-index+1
    cell(j,:)=strsplit(stringCell(j+index-1),' '); % split each line
end  
if(double(L(1))==0) % If the id of the first point (first line) of the cell is 0, increment the ids and the parents'id of every cell  to +1 
    cell(:,1)=cell(:,1)+1; % increment the ids of the points to +1 
    cell(:,7)=cell(:,7)+1; % increment the parents'ids to +1
end
cell=double(cell); % convert the cell into double
f=cell; %return the cell 
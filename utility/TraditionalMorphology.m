function f = TraditionalMorphology(Labels,XYZ_Coords,Parent_List)

% Create  directories where the DATA will be saved
if exist('TraditionalMorphology','file')==7 % If any the folder 'OUTPUTS' already exists, delete it
    rmdir('TraditionalMorphology','s');
end
mkdir('TraditionalMorphology'); % create a directory where the DATA about each cells will be saved

n = length(Labels);
X = XYZ_Coords(:,1);
Y = XYZ_Coords(:,2);
Z = XYZ_Coords(:,3);

Max_Euclidean_Distance=NeuronRadius(XYZ_Coords); 

A = sort(Parent_List);
Branching_Points = zeros(1,1);
Number_of_Stems = 0;
index = 1;
for i = 1:(length(A)-1)
    if A(i) == 1
        Number_of_Stems = Number_of_Stems + 1;
    end
    if (A(i+1)==A(i)) && (A(i) > 1)
        Branching_Points(index) = A(i);
        index = index + 1;
    end
end
Branching_Points = unique(Branching_Points);
Bifurcation_Number=length(Branching_Points);

Branching_Order = zeros(n,1);
Number_of_Branches = 0;
Path_Distance = zeros(n,1);
Total_Length = 0;
for i=1:n
    if (Parent_List(i) == 1) || (any(Branching_Points==Parent_List(i)))
        Branching_Order(i) = Branching_Order(Parent_List(i)) + 1;
        Number_of_Branches = Number_of_Branches + 1;
        Path_Distance(i) = Path_Distance(Parent_List(i))+((X(Parent_List(i))-X(i))^2+(Y(Parent_List(i))-Y(i))^2+(Z(Parent_List(i))-Z(i))^2)^0.5;
        Total_Length = Total_Length + ((X(Parent_List(i))-X(i))^2+(Y(Parent_List(i))-Y(i))^2+(Z(Parent_List(i))-Z(i))^2)^0.5;
    else if not(Parent_List(i) == -1) && not(Parent_List(i) == 0)
        Branching_Order(i) = Branching_Order(Parent_List(i));
        Path_Distance(i) = Path_Distance(Parent_List(i))+((X(Parent_List(i))-X(i))^2+(Y(Parent_List(i))-Y(i))^2+(Z(Parent_List(i))-Z(i))^2)^0.5;
        Total_Length = Total_Length + ((X(Parent_List(i))-X(i))^2+(Y(Parent_List(i))-Y(i))^2+(Z(Parent_List(i))-Z(i))^2)^0.5;
        end
    end
end
Maximum_Branch_Order = max(Branching_Order);
Maximum_Path_Distance = max(Path_Distance);

Tips = zeros(1,1);
index=1;
for i = 1:(length(A)-1)
    if not(Labels(i) == Parent_List(i+1))
        Tips(index) = Labels(i);
        index = index + 1;
    end
end
Number_of_Tips = length(Tips);

save('TraditionalMorphology/Number of Stems','Number_of_Stems');
save('TraditionalMorphology/Bifurcation Number','Bifurcation_Number');
save('TraditionalMorphology/Number of Branches','Number_of_Branches');
save('TraditionalMorphology/Number of Tips','Number_of_Tips');
save('TraditionalMorphology/Total Length','Total_Length');
save('TraditionalMorphology/Max Euclidean Distance','Max_Euclidean_Distance');
save('TraditionalMorphology/Maximum Path Distance to Root','Maximum_Path_Distance');
save('TraditionalMorphology/Maximum Branch Order','Maximum_Branch_Order');


clc;
%clear all;
%close all;

disp(' ## PROGRAM FOR CIRCLE DETECTION LEARNING AUTOMATA : PROJECT ## ');

%file directory
    
directory=char(pwd);
ImageDirectory = 'LA_Circle_Test_Images\';
ImageFiles = dir(ImageDirectory);
scale=0.5;
Edge_Pt_Num = 10000;

while Edge_Pt_Num>1200;   % Adjusting Image Size
    
    origIm=imread([ImageDirectory ImageFiles(10).name]);  % change the number in ImageFiles(3) from 3 to 10 for all images
    origIm=imresize(origIm,scale);
    origIm1=rgb2gray(origIm);
    BW = edge(origIm1,'canny');

    [rows, columns] = size(BW);
    Edge_Pt_Num = sum(BW(:));
    Edge_Pt_Data = zeros(Edge_Pt_Num,2);
    Radius_Range = [2, floor(sqrt((rows*rows)+(columns*columns))/2)]; % Rmin & Rmax for algorithm
    scale=scale-0.01;
    
end 

showIm=origIm;
imshow(showIm);
count=1;

% Creating Edge Point Data Matrix

for ii = 1:columns
    
    for i = 1:rows
        
        if BW(i,ii)==1 
            
            y=-i;
            x=ii;
       
            Edge_Pt_Data(count,1) = x;
            Edge_Pt_Data(count,2) = y;
            count=count+1;     
        end
    end
end

% Randomly Selecting 5 Percent Of Total Edge Points

Edge_Pt_Sel_Num=round((Edge_Pt_Num*5)/100);
Edge_Pt_Sel_Val=randperm(Edge_Pt_Num,Edge_Pt_Sel_Num);

% Forming Matrix For All Combinations Of Possible Candidate Circles (3 Pts)

Circle_Candidates_Val = nchoosek(Edge_Pt_Sel_Val,3);

% Computing Parameters For All The Circle Candidates (Center & Radius)

Circle_Candidates_Param=[];
Circle_Candidates_Param_final=[];

for num =1:size(Circle_Candidates_Val,1)
    
    % Co Ordinates Of 3 Non Colinear Points
    
    Pt_A = [Edge_Pt_Data(Circle_Candidates_Val(num,1),1), Edge_Pt_Data(Circle_Candidates_Val(num,1),2)];
    Pt_B = [Edge_Pt_Data(Circle_Candidates_Val(num,2),1), Edge_Pt_Data(Circle_Candidates_Val(num,2),2)];
    Pt_C = [Edge_Pt_Data(Circle_Candidates_Val(num,3),1), Edge_Pt_Data(Circle_Candidates_Val(num,3),2)];
    
    A = [(Pt_B(1)^2 + Pt_B(2)^2 - (Pt_A(1)^2 + Pt_A(2)^2)), 2*(Pt_B(2)-Pt_A(2)); (Pt_C(1)^2 + Pt_C(2)^2 - (Pt_A(1)^2 + Pt_A(2)^2)), 2*(Pt_C(2)-Pt_A(2))];
    B = [2*(Pt_B(1)-Pt_A(1)), (Pt_B(1)^2 + Pt_B(2)^2 - (Pt_A(1)^2 + Pt_A(2)^2)); 2*(Pt_C(1)-Pt_A(1)), (Pt_C(1)^2 + Pt_C(2)^2 - (Pt_A(1)^2 + Pt_A(2)^2))];
    
    Den = (4*((Pt_B(1)-Pt_A(1))*(Pt_C(2)-Pt_A(2))-(Pt_C(1)-Pt_A(1))*(Pt_B(2)-Pt_A(2))));
    X_Zero = floor(det(A)/Den);
    Y_Zero = floor(det(B)/Den);
    Radius = floor(sqrt(((X_Zero-Pt_A(1))^2) + ((Y_Zero-Pt_A(2))^2)));
    
    % Eliminating The Candidates Outside The Interested Radius Range
   
    if X_Zero>0.45*columns && X_Zero<0.55*columns && Y_Zero<(-0.45)*rows && Y_Zero>(-0.55)*rows && Radius>Radius_Range(1) && Radius<Radius_Range(2)
        
        Circle_Candidates_Param=[Circle_Candidates_Param; [X_Zero Y_Zero Radius]];
        
    end          
end

Circle_Candidates_Param = unique(Circle_Candidates_Param,'rows');

% Eliminating The Candidates Which Are Similar To Each Other 

num=1;

while num <= size(Circle_Candidates_Param,1)
    
    X_Zero_Range = ((Circle_Candidates_Param(num,1)-1):(Circle_Candidates_Param(num,1)+1));
    Y_Zero_Range = ((Circle_Candidates_Param(num,2)-1):(Circle_Candidates_Param(num,2)+1));
    Radius_Range1 = ((Circle_Candidates_Param(num,3)-1):(Circle_Candidates_Param(num,3)+1));
    Repeat_Param_Vector = combvec(X_Zero_Range,Y_Zero_Range,Radius_Range1);
    Repeat_Param_Vector = Repeat_Param_Vector';
    %Repeat_Param_Vector = removerows(Repeat_Param_Vector,'ind',[round(size(Repeat_Param_Vector,1)/2)]);    
    X_Zero=Circle_Candidates_Param(num,1);
    Y_Zero=Circle_Candidates_Param(num,2);
    Radius=Circle_Candidates_Param(num,3);
    
    Circle_Candidates_Param(any(ismember(Circle_Candidates_Param,Repeat_Param_Vector,'rows'),2),:) = []; 
    Circle_Candidates_Param_final = [Circle_Candidates_Param_final; [X_Zero Y_Zero Radius]];
    
    num=num+1;
end


Remaining_Actions = size(Circle_Candidates_Param_final,1);
K_max = round(Remaining_Actions/2);

% Learning Automata Algorithm

Inital_Probability = (1/Remaining_Actions);
Prob_Update_Matrix = Inital_Probability*(ones(Remaining_Actions,1));
Prob_Action_Matrix = Prob_Update_Matrix;
%Beta = zeros(Remaining_Actions,1);
Learning_Rate = 0.005;


for i=1:K_max
    
     disp(i);
    
     % 1. Selectiion of actions 'Av' for the iteration 

     z=rand; % Selecting Pseudo Random Number Between 0-1
     Prob_Update_Step = [];
     count=1;
     n=randperm(Remaining_Actions)';

     while sum(Prob_Update_Step)<z
    
          Prob_Update_Step(count,1)=Prob_Update_Matrix(n(count,1),1);
          count=count+1;
          
     end
     
     if count>=3
         
         n = n(1:(count-2),1);
         Prob_Update_Step = Prob_Update_Step(1:(count-2),1);
     
     end

     % 2. Calculation of reinforcement signal for the actions 'Av'
     
     Beta_Update_Vector=zeros(size(n,1),2);

     for ii = 1:size(n,1)
     
         [x1, y1] = getmidpointcircle(Circle_Candidates_Param_final(n(ii,1),1), Circle_Candidates_Param_final(n(ii,1),2), Circle_Candidates_Param_final(n(ii,1),3)); 
         Error_Beta_Cal=zeros(size(x1));
    
    
         for m = 1:size(x1,1)
        
             if ismember([x1(m,1), y1(m,1)],Edge_Pt_Data,'rows') == true
            
                  Error_Beta_Cal(m,1) = 1;
             
             end
        
         end
    
         Beta_Cal = (sum(Error_Beta_Cal))/size(x1,1);
         Beta_Update_Vector(ii,1) = n(ii,1);
         Beta_Update_Vector(ii,2) = Beta_Cal;
         %Beta(n(ii,1),1) = Beta_Cal;
    
         %Prob_Update_Matrix(n(ii,1),1) = Prob_Update_Matrix(n(ii,1),1) + (Learning_Rate*Beta(n(ii,1),1)*(1-Prob_Update_Matrix(n(ii,1),1)));
    
     end
     
     % Probability Update Law
     
     [row_max] = find(Beta_Update_Vector(:,2)==max(Beta_Update_Vector(:,2)));
     Action_Num_Max = Beta_Update_Vector(row_max,1);
     Beta_Max = Beta_Update_Vector(row_max,2);
     
     % A. Rewarding The Action With Highest Reinforcement Signal
     
     %for ii = 1:size(row_max,1)
         
          %Prob_Update_Matrix(Action_Num_Max(ii,1),1) = Prob_Update_Matrix(Action_Num_Max(ii,1),1) + (Learning_Rate*Beta_Max(ii,1)*(1-Prob_Update_Matrix(Action_Num_Max(ii,1),1)));
          %Beta_Update_Vector = removerows(Beta_Update_Vector,'ind',[row_max(ii,1)]);
     
      Prob_Update_Matrix(Action_Num_Max(1,1),1) = Prob_Update_Matrix(Action_Num_Max(1,1),1) + (Learning_Rate*Beta_Max(1,1)*(1-Prob_Update_Matrix(Action_Num_Max(1,1),1)));
      Beta_Update_Vector = removerows(Beta_Update_Vector,'ind',[row_max]);
         
          
     %end
     
     % B. Penalizing Rest Of The Actions
     
     for ii = 1:(size(n,1)-size(row_max,1))
         
         Action_Num = Beta_Update_Vector(ii,1);
         Beta_Val = Beta_Update_Vector(ii,2);
         
         Prob_Update_Matrix(Action_Num,1) = Prob_Update_Matrix(Action_Num,1) - (Learning_Rate*Beta_Val*Prob_Update_Matrix(Action_Num,1));
         
     end
  
     Prob_Action_Matrix = [Prob_Action_Matrix, Prob_Update_Matrix];

end

surf(Prob_Action_Matrix);
title(' Evolution of Probability Distribution Curve ');
xlabel(' Iteration Number ');
ylabel(' Action Number ');
zlabel(' Probability Value ');
figure;

% plotting estimated circle

Final_Action_Prob=max(Prob_Action_Matrix(:,size(Prob_Action_Matrix,2)));
[action_row_max] = find(Prob_Action_Matrix(:,size(Prob_Action_Matrix,2))==Final_Action_Prob);

[xf, yf] = getmidpointcircle(Circle_Candidates_Param_final(action_row_max,1),(rows+Circle_Candidates_Param_final(action_row_max,2)),Circle_Candidates_Param_final(action_row_max,3));


showIm=origIm;
imshow(showIm);
hold on;
plot(xf,yf,'r');
title(' Predicted Circle (Red) With Most Confidence ');

% For analysis:

[sortedX,sortingIndices] = sort(Prob_Action_Matrix(:,size(Prob_Action_Matrix,2)),'descend');

figure;
pks = findpeaks(Prob_Action_Matrix(:,size(Prob_Action_Matrix,2)'));
plot(pks);
title(' Peaks of Probability At Final Iteration ');
xlabel(' Peak Number ');
ylabel(' Probability Value ');

max_val=max(sortedX);


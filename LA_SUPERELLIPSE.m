clc;
%clear all;
%close all;

disp(' ## PROGRAM FOR TREE RINGS DETECTION LEARNING AUTOMATA : PROJECT ## ');

%file directory
    
directory=char(pwd);
ImageDirectory = 'LA_Superellipse_Test_Images\';
ImageFiles = dir(ImageDirectory);
scale=0.5;
Edge_Pt_Num = 10000;

while Edge_Pt_Num>2500;   % Adjusting Image Size
    
    origIm=imread([ImageDirectory ImageFiles(13).name]);    % change the number in ImageFiles(3) from 3 to 13 for all images
    origIm=imresize(origIm,scale);                         % for ImageFiles(3),ImageFiles(4),ImageFiles(5) take N=2 and K=0.5 below
    origIm1=rgb2gray(origIm);                              % for rest of the images take N=1.9 and K=0.95 below
    BW = edge(origIm1,'canny');
    disp(Edge_Pt_Num);
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
Edge_Pt_Data1=Edge_Pt_Data;
Edge_Pt_Data1=removerows(Edge_Pt_Data1,'ind',Edge_Pt_Sel_Val');
Edge_Pt_Num1=size(Edge_Pt_Data1,1);
Edge_Pt_Sel_Val1=randperm(Edge_Pt_Num1,Edge_Pt_Sel_Num);

% Forming Matrix For All Combinations Of Possible Superellipse Candidates (6 Pts)

Circle_Candidates_Val1 = nchoosek(Edge_Pt_Sel_Val,3);
Circle_Candidates_Val2 = nchoosek(Edge_Pt_Sel_Val1,3);

Circle_Candidates_Val = [Circle_Candidates_Val1, Circle_Candidates_Val2];


% Computing Parameters For All The Circle Candidates (Center & Radius)

Circle_Candidates_Param=[];
Circle_Candidates_Param_final=[];
warning('off');


for num =1:size(Circle_Candidates_Val,1)
    
    points=Circle_Candidates_Val(num,:);
    
    if size(points,2)==size(unique(points),2) 
        
        Xp=[Edge_Pt_Data(points(1,1),1);Edge_Pt_Data(points(1,2),1);Edge_Pt_Data(points(1,3),1);Edge_Pt_Data1(points(1,4),1);Edge_Pt_Data1(points(1,5),1);Edge_Pt_Data(points(1,6),1)];
        Yp=[Edge_Pt_Data(points(1,1),2);Edge_Pt_Data(points(1,2),2);Edge_Pt_Data(points(1,3),2);Edge_Pt_Data1(points(1,4),2);Edge_Pt_Data1(points(1,5),2);Edge_Pt_Data(points(1,6),2)];
        E=fit_ellipse(Xp,Yp);
    
         X_Zero=E.X0_in;
         Y_Zero=E.Y0_in;
         Radius=round(E.long_axis/2);
    
        % Eliminating The Candidates Outside The Interested Radius Range
   
        if sum(X_Zero)~=0 || sum(Y_Zero)~=0 || sum(Radius)~=0
    
            if X_Zero>0.46*columns && X_Zero<0.54*columns && Y_Zero<(-0.46)*rows && Y_Zero>(-0.54)*rows && Radius>Radius_Range(1) && Radius<Radius_Range(2)
            
                Circle_Candidates_Param=[Circle_Candidates_Param; [X_Zero Y_Zero Radius]];
        
            end
        end
    end
end

Circle_Candidates_Param = unique(Circle_Candidates_Param,'rows');
Circle_Candidates_Param=floor(Circle_Candidates_Param);

% Eliminating The Candidates Which Are Similar To Each Other 

num=1;

while num <= size(Circle_Candidates_Param,1)
    
    X_Zero_Range = ((Circle_Candidates_Param(num,1)-1):(Circle_Candidates_Param(num,1)+1));
    Y_Zero_Range = ((Circle_Candidates_Param(num,2)-1):(Circle_Candidates_Param(num,2)+1));
    Radius_Range1 = ((Circle_Candidates_Param(num,3)-2):(Circle_Candidates_Param(num,3)+2));
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
Learning_Rate = 0.0005;

% Ellipse properties:

N=1.9;
K=0.95;



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
     
         [x1, y1] = getsuperellipse(Circle_Candidates_Param_final(n(ii,1),1), Circle_Candidates_Param_final(n(ii,1),2), Circle_Candidates_Param_final(n(ii,1),3),N,K); 
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
     
     % A. Rewarding The Actions With Highest Reinforcement Signal
     
     for ii = 1:size(row_max,1)
         
          Prob_Update_Matrix(Action_Num_Max(ii,1),1) = Prob_Update_Matrix(Action_Num_Max(ii,1),1) + (Learning_Rate*Beta_Max(ii,1)*(1-Prob_Update_Matrix(Action_Num_Max(ii,1),1)));
          %Beta_Update_Vector = removerows(Beta_Update_Vector,'ind',[row_max(ii,1)]);
     
      %Prob_Update_Matrix(Action_Num_Max(1,1),1) = Prob_Update_Matrix(Action_Num_Max(1,1),1) + (Learning_Rate*Beta_Max(1,1)*(1-Prob_Update_Matrix(Action_Num_Max(1,1),1)));
      %Beta_Update_Vector = removerows(Beta_Update_Vector,'ind',[row_max]);
         
          
     end
     
     Beta_Update_Vector = removerows(Beta_Update_Vector,'ind',[row_max]);
     
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

showIm=origIm;
imshow(showIm);
hold on;


for i=1:size(action_row_max,1)
    
   [xf, yf] = getsuperellipse(Circle_Candidates_Param_final(action_row_max(i,1),1),(rows+Circle_Candidates_Param_final(action_row_max(i,1),2)),Circle_Candidates_Param_final(action_row_max(i,1),3),N,K);
   plot(xf,yf,'r');
   title(' Predicted Tree Ring(Red) With Most Confidence ');
   hold on;

end

% For analysis:

[sortedX,sortingIndices] = sort(Prob_Action_Matrix(:,size(Prob_Action_Matrix,2)),'descend');

max_val=max(sortedX);
Tree_Ring_Num=0;

for val=1:size(sortedX,1);
    
    peak_val=sortedX(val,1);
    
    if peak_val > max_val*0.50; % Peaks below 50% of max value are rejected
        
        Tree_Ring_Num=Tree_Ring_Num+1;
     
    end
        
end


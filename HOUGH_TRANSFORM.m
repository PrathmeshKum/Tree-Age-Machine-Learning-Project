clc;
clear all;
close all;

disp(' ## PROGRAM FOR CIRCULAR HOUGH TRANSFORM : PROJECT ## ');

%file directory
    
directory=char(pwd);
ImageDirectory = 'Hough_Test_Images\';
ImageFiles = dir(ImageDirectory);

origIm=imread([ImageDirectory ImageFiles(7).name]);  % change the number in ImageFiles(3) from 3 to 7 for all images

%origIm=imresize(origIm,0.5);  % scale image appropriately if neccessary
origIm1=rgb2gray(origIm);      % converting into grayscale
BW = edge(origIm1,'canny');    % getting edge points
showIm=origIm;
imshow(showIm);


[rows, columns]=size(BW);

% Creating parameter space

RADIUS=(25:floor(sqrt((rows*rows)+(columns*columns))/2));    % Select appropriate range of radius to search
THETA=(0:90);
A=(1:floor(3*columns));
B=(1:floor(3*rows));
PARAMETER_SPACE=zeros(size(A,2),size(B,2),size(RADIUS,2));   

for i = 1:rows
    
    for ii = 1:columns
        
        if BW(i,ii)==1 
            
            x=ii;
            y=-i;
            
            for r = RADIUS
                
            
                for t = 1:(size(THETA,2))
                
                    a(1,t)=(x-r*cosd(-(t-1)))+(columns);   % Shifting the values by factor of 1
                    b(1,t)=-(y+r*sind(-(t-1)))+(rows);     % Parametric equations taken as mentioned in report
                
                    % updating votes in parameter space
                
                
                    PARAMETER_SPACE(floor(a(1,t)),floor(b(1,t)),(r-(RADIUS(1)-1)))=PARAMETER_SPACE(floor(a(1,t)),floor(b(1,t)),(r-(RADIUS(1)-1)))+1;
                
                end 
            end
         end
        
    end
 
end

final_point_count=max(PARAMETER_SPACE(:));  % Determining the detected circle
final_radius=[];
 
for c= (RADIUS-(RADIUS(1)-1))
    
    % Obtaining parameter values for detected circle
    
    d=PARAMETER_SPACE(:,:,c);
    point_count=max(d(:));
    
    if point_count==final_point_count;
        
        final_radius=[final_radius; c];
        [aa,bb]=find(d==final_point_count);
        aa=aa-(columns);
        bb=-(bb-(rows));
        
        break
        
    else
        continue
    end
    
end

final_radius=max(final_radius)+(RADIUS(1)-1);
aa=max(aa);
bb=max(bb);
figure;
mesh(d);    % Plotting the parameter surface for determined circle
title('Parameter Space Plot for Radius of Detected Circle');
xlabel(' X coordinate ');
ylabel(' Y coordinate ');
zlabel(' Number of votes ');

for t = 1:(size(THETA,2)*4)
    
     xunit(1,t) = final_radius*cosd(t) + aa;
     yunit (1,t)= final_radius*sind(t) + bb;
   
end

figure;
h = plot(xunit,yunit,'r'); % plotting the detected circle
title('Detected Circle Plot');
axis([0 columns -rows 0]);

% Equation for plotting shape of tree ring

% Superellipse

theta=(0:360);
radius = [];
b=(2:4:50);
n=1.6;
k=0.96;
%a=50;
x1=10;
y1=5;

for a = b
    
     for i = 1:size(theta,2)
    
          radius(1,i) = a*(((abs(cosd(theta(i)))^n) + (abs(sind(theta(i))/k)^n))^(-1/n));
          x(1,i) = x1+(radius(1,i)*cosd(theta(i)));
          y(1,i) = y1+(radius(1,i)*sind(theta(i)));
    
     end

     plot(x,y);
     hold on;
end

plot(x1,y1,'*');
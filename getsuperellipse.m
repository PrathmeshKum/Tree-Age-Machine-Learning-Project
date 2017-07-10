function [x1, y1] = getsuperellipse(xcen,ycen,a,N,K)

 theta=(0:360);
 
 radius = a.*(((abs(cosd(theta)).^N) + (abs(sind(theta)/K).^N)).^(-1/N));
 x1 = xcen+(radius.*cosd(theta));
 y1 = ycen+(radius.*sind(theta));
 
end
%--------------------------------------------------------------------------
%     GP
%--------------------------------------------------------------------------
function gp=GP(t,colpoint,sourceregion,prop,freq)
%xe,ye,nx,ny,i)     
xsp =sourceregion.x(2)-sourceregion.x(1);
ysp =sourceregion.y(2)-sourceregion.y(1);
xs  =sourceregion.x(1)+t*xsp;
ys  =sourceregion.y(1)+t*ysp;
r   =norm([colpoint.x-xs,colpoint.y-ys]);
difs=sqrt(xsp^2+ysp^2);
k=freq.omega/prop.c;
gp=-.25*i*besselh(0,2,k*r)*difs;
return
end            

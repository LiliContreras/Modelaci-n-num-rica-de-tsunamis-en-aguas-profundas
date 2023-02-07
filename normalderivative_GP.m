%--------------------------------------------------------------------------
%     FUNCIÓN PARCIAL DE GP CON RESPECTO A LA NORMAL.
%--------------------------------------------------------------------------
function dgp=normalderivative_GP(t,colpoint,sourceregion,vnormal,prop,freq)
    
xsp =sourceregion.x(2)-sourceregion.x(1);
ysp =sourceregion.y(2)-sourceregion.y(1);
xs  =sourceregion.x(1)+t*xsp;
ys  =sourceregion.y(1)+t*ysp;
r   =norm([colpoint.x-xs,colpoint.y-ys]);
gx  =(colpoint.x/r)-xs/r;
gy  =(colpoint.y/r)-ys/r;
difs=sqrt(xsp^2+ysp^2);
k=freq.omega/prop.c;
dgp=.25*i*k*besselh(1,2,k*r)*((gx*vnormal.n(1))+(gy*vnormal.n(2)))*difs;
return
end            

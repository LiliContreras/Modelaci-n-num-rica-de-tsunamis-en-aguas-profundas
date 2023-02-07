%--------------------------------------------------------------------------
%     Integración gaussiana usando 10 puntos
%--------------------------------------------------------------------------
function inteicis=inte_ic_is(colpoint,vnormal,sourceregion,prop,freq)
%Coeficientes
w=[.2955242247,.2692667193,.2190863625,.1494513491,.0666713443];
x=[.1488743389,.4333953941,.6794095682,.8650633666,.9739065285];

xM=0.5; %Puntos medios de segmentos
xR=0.5;
inte=0;
for J=1:5
  dx=xR*x(J) ;
    inte=inte+w(J)*(normalderivative_GP(xM+dx,colpoint,sourceregion,vnormal,prop,freq)+...
                    normalderivative_GP(xM-dx,colpoint,sourceregion,vnormal,prop,freq));
end
inteicis=xR*inte;
return
end


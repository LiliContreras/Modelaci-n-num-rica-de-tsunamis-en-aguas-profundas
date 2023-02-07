%
% Cálculo de refinamiento
%
function [xN,yN,nx,ny,s]=refine_boundarymesh(x,y,ds,condiciones)%meshfactor)
% x and y rectified to have the same length
% compute the acumulated distance s

%Se construye refinamiento dependiendo del método en que se establecen condiciones
if condiciones==0 %Método de imágenes
	s(1)=0;
	for iP=2:length(x)
		%Se calcula distancia acumulada de S (distancia de superficie)
    		s(iP)=s(iP-1)+norm([x(iP) y(iP)]-[x(iP-1) y(iP-1)]);
	end
	%Segmentos en los que se divide la superficie
	nP=1+round(s(end)/ds);
	%Construcción de coordenadas X,Y
	sV=linspace(x(2),s(iP),nP-1);
	xN=interp1(s,x,sV);
	xNi=[0];
	xN=[xNi,xN];
	yN=interp1(s,y,sV);
	yNi=[y(1)];
	yN=[yNi,yN];
	%Vectir auxiliar para cálculo de vector normal
	k=[0  0 1];
	%Cálculo de vector normal a segmento de superficie
	for iP=2:length(xN) 
    		vn=cross([xN(iP) yN(iP) 0 ]-[xN(iP-1) yN(iP-1) 0],k);
    		vn=-vn/norm(vn);
    		nx(iP-1)=vn(1);
    		ny(iP-1)=vn(2);
	end
else
%Se construye modelo con condiciones en superficie
	s(1)=0;
	for iP=2:length(x)
    		s(iP)=s(iP-1)+norm([x(iP) y(iP)]-[x(iP-1) y(iP-1)]);
	end

	nP=1+round(s(end)/ds);


	sV=linspace(x(2),s(iP),nP-1);
	xN=interp1(s,x,sV);
	for i=1:length(xN)
    		if xN(i)<max(x)
        		num=i-1; %contador auxiliar
    		end
	end
	xNi=[0];
	xNs=fliplr(xN(1:num));
	xN=[xNi,xN,xNs,0];
	yN=interp1(s,y,sV);
	yNi=[y(1)];
	yNs=zeros(1,length(xNs));
	yN=[yNi,yN,yNs,0];
	k=[0  0 1];

	for iP=2:length(xN) 
    		vn=cross([xN(iP) yN(iP) 0 ]-[xN(iP-1) yN(iP-1) 0],k);
    		vn=-vn/norm(vn);
    		nx(iP-1)=vn(1);
    		ny(iP-1)=vn(2);   
	end
end

end

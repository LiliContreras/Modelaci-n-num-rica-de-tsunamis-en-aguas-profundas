%
% Se construye imágen del modelo
%

function boundary=build_image(x,y,nx,ny)


nP=2*length(x)-1;

for i=1:length(x)-1
    xN(i)=x(i); %Original
    xN(nP-i+1)=x(i); %Imagen
    yN(i)=y(i);
    yN(nP-i+1)=-y(i);
    %Vectores normales de modelo original y de la imagen
    nxN(i)=nx(i);
    nyN(i)=ny(i);
    nxN(nP-i)= nx(i);
    nyN(nP-i)=-ny(i);
end
xN(length(x))=x(end);
yN(length(x))=y(end);
%Puntos centrales entre segmentos
xCN=.5*(xN(1:end-1)+xN(2:end));
yCN=.5*(yN(1:end-1)+yN(2:end));  

%Se cargan resultados en boundary
boundary.xN  = xN;
boundary.yN  = yN;
boundary.nxN = nxN;
boundary.nyN = nyN;
boundary.xCN = xCN;
boundary.yCN = yCN;



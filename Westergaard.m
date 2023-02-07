%Solución analítica de presiones sobre contacto lateral, Westergaard
function [PPMAXw YHw]=Westergaard(f)
%========================================================================
%Propiedades
rho=1000; %Densidad
f=f;  %Frecuencia
W=2*pi()*f; %Frecuencia angular
U=1; %Desplazamiento unitario
C=1500; %Velocidad
K=W/C; %Número de onda
h=5000; %Profundidad máxima
N=20;  %Número de elementos

y=linspace(0,5000,11); %Vector de profundidad
x=ones(1,length(y))*1; %Coordenadas x
P=zeros(1,length(y)); %Vector de Presiones
for i=1:length(y) %Se construye función analítica
    for n=1:N
        a=(-1)^(n+1);
        b=(2*n-1);
        lambda=b*pi()/(2*h);
        sn=sqrt(K^2-lambda^2);
        Q=-(4*a*rho*W^2*U)/(pi()*b*sn);
        P(i)=P(i)+Q*cos(lambda*y(i))*exp(-i*sn*x(i)); %Cálculo de presión para posición i
    end
end

%Cálculo relación Presión máxima y profundidad
PPMAXw=abs(P/max(abs(P)));
PPMAXw=fliplr(PPMAXw);
YHw=y/h;
YHw=fliplr(YHw);

end

%======================================================================
%	  Modelación numérica de tsunamis en aguas profundas
% Aplicación de la teoría de ecuaciones integrales para un modelo de
% aguas profundas.

% Lilibeth Zaira Contreras Alvarado
% Leonardo Ramírez Guzmán
%======================================================================

%=========================PARÁMETROS DEL MÓDELO=========================
% Definición de frecuencias
frfinal=.125;
frinic=.01;
%vector de frecuencias
f=linspace(frinic,frfinal,64); 
df=f(2)-f(1);
%Número de puntos por segmento
npplo=4;

%velocidad 
c=1500;
prop.c=c;
%densidad del agua
rho=1000;
prop.rho=rho;
%número de puntos de observación
nest=100;

%Batimetría Modelo 1:
%x=[ 0 20 160 160.1]*1e3; %50 metros, 150000
%y=[-5 -5 -5  0]*1e3; %30metros

%Batimetria Modelo 2:
x=[0 20 300 700 700.1]*1e3;
y=[-5 -5 -5 0 0]*1e3;

%Puntos de observación en superficie
%Modelo 1
%observer.xO=linspace(0,260,nest)*1e3; %metros
%observer.yO=0*observer.xO;

%Modelo 2
observer.xO=linspace(0,1000,nest)*1e3; %metros
observer.yO=0*observer.xO;

%Se crea malla de batimetría
%	Condiciones:
% Imágenes=0; PresionCero=1; Gravedad=2;
C=2;

[xN,yN,nx,ny,S]=refine_boundarymesh(x,y,1e4,C);

if C~=0
	xC=.5*(xN(1:end-1)+xN(2:end));
	yC=.5*(yN(1:end-1)+yN(2:end));
	boundary.xN  = xN;
	boundary.yN  = yN;
	boundary.nxN = nx;
	boundary.nyN = ny;
	boundary.xCN = xC;
	boundary.yCN = yC;
else
	boundary=build_image(xN,yN,nx,ny);
end

%Grafica de segmentos
figure;
plot(x,y)
hold on
plot(boundary.xN,boundary.yN,'r.','LineStyle','none')
legend('Contorno del modelo','Separación de segmentos','FontSize', 13)
ylabel('Profundidad [m]','FontSize', 14)
xlabel('Distancia x [m]','FontSize', 14)

%variable auxiliar, da índice de donde comienza superficie del modelo
for i=1:length(yN)
    if yN(i)<0
        aux=i+1;
    end
end

%Puntos de observación en contacto lateral
yO=[0:(yN(aux)-yN(aux-1))/2:5000];
observerl.xO=ones(1,length(yO))*16e4;
observerl.yO=-yO;
%=======================================================================


%================================FUENTE=================================
%Fuente 1: Desplazamiento lateral F=1;
theta=0;
disprig=[cosd(theta) sind(theta)];
disprig=disprig/norm(disprig);

%Fuente 2: Desplazamiento vertical unitario con retraso aleatorio F=2;
disp=[0,0,0,0,0,0,0,ones(1,length(f)-7)];
uza=zeros(round(aux/2)-1,length(f));
uzb=zeros(round(aux/2),length(f));
disp=[uza;disp;uzb];

%Fuente 3: Desplazamiento vertical, Okada F=3;
%Se carga desplazamiento vertical Okada
DVertical=importdata('DespVertical.txt');
a=DVertical.data(:,2)';
uz=[a,zeros(1,aux-length(a))];
%Construcción de función rampa
%Tiempo de corte
t0=30;%s
%retraso de tiempo
rt=20; %s
%Vector de tiempo
t=[0:1:length(f)-1-rt];
%Variables auxiliares para almacenar resultados
r=ones(length(uz),length(t));
c=zeros(1,length(uz));
d=zeros(1,length(uz));
rampa=zeros(length(uz),length(d)+1);
for j=1:length(uz)
    for i=1:length(t)
        if t(i)<=0 %se pone 0 antes del tiempo de recepción
            r(j,i)=0;
            c(j)=c(j)+1;
        elseif t(i)>=t0 %se asigna el valor mayor después de tiempo de corte
            r(j,i)=uz(1,j);
        else
            d(j)=d(j)+1; %contador auxiliar
        end
    end
end
%Se llenan valores de la rampa
for i=1:length(uz)
    for j=1:d(i)
        rampa=linspace(0,uz(1,i),d(i)+1);
        r(i,c(i)+j)=rampa(j+1);
    end
end
R=[zeros(aux,rt),r]; %Se agrega retraso de tiempo
for i=1:length(uz)
	%Se transforma al dominio de frecuencias
    Fuz(:,i)=fftshift(fft(R(i,:)));    
end

%Se define la fuente a usar:
F=3;
%========================================================================


%======================CÁLCULO DE INCOGNITA XI==========================
%Se obtiene incognita para cada frecuencia
for i=1:length(f)
    %Longitud de onda
    lambda=prop.c/f(i);
    dS=lambda/npplo; 
    %Se crean segmentos para cada frecuencia
    [xNf,yNf,nxf,nyf,Sf]=refine_boundarymesh(x,y,min(dS,1e4),C);
    
	if C~=0
		xCf=.5*(xNf(1:end-1)+xNf(2:end));
		yCf=.5*(yNf(1:end-1)+yNf(2:end));
		boundary.xN  = xNf;
		boundary.yN  = yNf;
		boundary.nxN = nxf;
		boundary.nyN = nyf;
		boundary.xCN = xCf;
		boundary.yCN = yCf;
	else
		%Se construye imagen del modelo
		boundary=build_image(xNf,yNf,nxf,nyf);
	end
    
    %Frecuencia y frecuencia angular
    freq.f=f(i);
    freq.omega=2*pi*f(i);
    %FUENTE
    if F==1
	desp=disprig;
    elseif F==2
	desp=disp(:,i)';
    else
	desp=Fuz(i,:);
    end
    %Se construye sistema de ecuaciones y obtiene solución xi
    solution(i).Xi=build_solve_ie(boundary,prop,freq,desp,aux,C);    
    display('freq')
    f(i)
end
%========================================================================

if F==1
    pl="1-Cálculo de presión en pared lateral \n2-Cálculo de desplazamiento en superficie\n";
    opc=input(pl)
    if opc==1
        obs=observerl;
    else
        obs=observer;
    end
else
    opc=0;
    obs=observer;
end

%=======================CÁLCULO DESPLAZAMIENTO===========================
for i=1:length(f)
    lambda=prop.c/f(i);
    dS=lambda/npplo;    
    [xNf,yNf,nxf,nyf,Sf]=refine_boundarymesh(x,y,min(dS,1e4),C);  
    
    if C~=0
		xCf=.5*(xNf(1:end-1)+xNf(2:end));
		yCf=.5*(yNf(1:end-1)+yNf(2:end));
		boundary.xN  = xNf;
		boundary.yN  = yNf;
		boundary.nxN = nxf;
		boundary.nyN = nyf;
		boundary.xCN = xCf;
		boundary.yCN = yCf;
	else
		%Se construye imagen del modelo
		boundary=build_image(xNf,yNf,nxf,nyf);
	end

    freq.f=f(i);
    freq.omega=2*pi*f(i);
    freq.i=i;

    %Se obtiene desplazamiento y presiones en puntos de observación
    [un(i).p un(i).dispn]=compute_unormal(boundary,prop,freq,solution,obs,C);  
end
%========================================================================

%===================TRANSFORMANDO A DOMINIO DEL TIEMPO===================
for i=1:length(f)
   for j=1:length(obs.xO)
      Ad(i,j)=un(i).dispn(j); 
      Ap(i,j)=un(i).p(j); 
   end
end

figure;
surf(obs.xO,f,abs(Ad));
view(2)
shading interp
xlabel('x')
ylabel('frecuencia (Hz)')
title('Función de transferencia')
colorbar


UN=Ad*0;
for i=1:length(obs.xO)
    UN(:,i)=real(ifft(fftshift(Ad(:,i))));
end

%Cálculo de Presión en pared lateral 
if opc==1
    %Distribución de presión en contacto lateral
    %Para frecuencia 1
    p1=Ap(1,:);
    pmax=max(abs(p1));
    ppmax=abs(p1/pmax);

    h=5000;
    yp=yO/h;
    yp=fliplr(yp);

    [PPMAXw YHw]=Westergaard(0.01);

    figure
    plot(ppmax,yp,LineWidth=2)
    hold on
    plot(PPMAXw,YHw,'*',LineWidth=2)
    legend('Cálculo de presiones', 'Solución analítica')
else
    %figure
    %for i=1:length(f)
    %    plot(observer.xO,UN(i,:))
    %    ylim([-5 5])
    %    pause(0.3)
    %end
    figure
    for i=1:100
        xo=ones(1,64)*i;
        t=linspace(0,100,64);
        hold on
        plot3(xo,t,UN(:,i),'b')
        zlim([-1 1])
    end
end

%{
figure
filename = 'SimulacionAlturaOla_DespVertical_ConGravedad.gif';
for n = 1:length(f)
      plot(observer.xO,UN(n,:))
      ylim([-2 2])
      drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if n == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end
end
%}


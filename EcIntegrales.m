%======================================================================
%	  Modelación numérica de tsunamis en aguas profundas
% Aplicación de la teoría de ecuaciones integrales para un modelo de
% aguas profundas.

% Lilibeth Zaira Contreras Alvarado
% Leonardo Ramírez Guzmán
%======================================================================
clear all
%close all
clc

%MODELO 1=1     MODELO 2=2
M=1;

%%
%=========================PARÁMETROS DEL MÓDELO=========================
% Definición de frecuencias
frfinal=0.5;%
frinic=0.01; %

if M==1
    tt=600;%s
else
    tt=900; %s
end
w0=2*pi/tt;
df=w0/(2*pi);          % Diferencial de frecuencia
f=[frinic:df:frfinal]; % Vector de frecuencias

%Número de puntos por segmento
npplo=20;

%velocidad 
c=1500; %m/s
prop.c=c;
%densidad del agua
rho=1000;
prop.rho=rho;
%número de puntos de observación
nest=100;

if M==1
    %Batimetría Modelo 1:
    x=[ 0 20 160 160.1]*1e3; % m
    y=[-5 -5 -5  0]*1e3;     % m
else
    %Batimetria Modelo 2:
    x=[0 20 300 700 700.1]*1e3; % m
    y=[-5 -5 -5 0 0]*1e3;       % m
end

%Se crea malla de batimetría
%	Condiciones:
% Imágenes=0; PresionCero=1; Gravedad=2;
C=0;     % Condición
DS=1e3;  % Resolución
[xN,yN,nx,ny,S]=refine_boundarymesh(x,y,DS,C);

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

%Grafica de segmentos de batimetría
figure;
plot(x/1000,y,LineWidth=2) % km 
hold on
plot(boundary.xN/1000,boundary.yN,'r.','LineStyle','none')
legend('Contorno del modelo','Separación de segmentos','FontSize', 13)
ylabel('Profundidad [m]','FontSize', 14)
xlabel('Distancia x [km]','FontSize', 14)

%Variable auxiliar, da índice de donde comienza superficie del modelo
for i=1:length(boundary.yCN)
    if boundary.yCN(i)<0
        aux=i+1; %Índice donde comienza superficie
    end
end
for i=1:length(boundary.yCN)
    if boundary.yCN(i)==-max(abs(boundary.yCN))
        aux2=i+1; %índice del fondo del modelo
    end
end

%Puntos de observación en superficie
if M==1
    %Modelo 1
    observer.xO=linspace(0,160.1,nest)*1e3;    % m
    observer.yO=zeros(1,length(observer.xO));  % m
else
    %Modelo 2
    observer.xO=linspace(0,700.1,nest)*1e3; % m
    observer.yO=0*observer.xO;              % m
end

%Puntos de observación para análisis de presiones
% en contacto lateral
yO=[0:(yN(aux)-yN(aux-1))/2:5000];
observerl.xO=ones(1,length(yO))*16e4;
observerl.yO=-yO;

%Puntos de observación para análisis de propagación
%Malla de observación
xM=linspace(0,xN(aux2),aux2);
yM=linspace(0,yN(aux2),10);
[xx yy]=meshgrid(xM,yM);
%=======================================================================

%%
%================================FUENTE=================================
%Aquí se muestra como se crea la fuente de manera general
%Pero la fuente a usar, se crea dentro del ciclo
%Fuente 1: Desplazamiento lateral F=1;
theta=0;
disprig=[cosd(theta) sind(theta)];
disprig=disprig/norm(disprig);

%Fuente 2: Desplazamiento vertical unitario F=2;
disp=ones(1,length(f));
uza=zeros(round(aux2/2)-1,length(f));
uzb=zeros(round(aux2/2),length(f));
uzl=zeros(aux-aux2,length(f));
disp=[uza;disp;uzb;uzl];

%Fuente 3: Desplazamiento vertical, Okada F=3;
%FUNCION RAMPA
DVertical=importdata('DespVertical.txt');
a=DVertical.data(:,2)';
uz=[a,zeros(1,aux-length(a))];
rt=50;  % 50 s de retraso de tiempo
t0=150; % 150 Tiempo de corte
trm=linspace(0,tt,length(f)*2); % Vector de tiempo
r=ones(1,length(f)*2);          % Se inicializa rampa en cero
c=0;
d=0;
for i=1:length(trm)
        if trm(i)<=rt %se pone 0 antes del tiempo de recepción
            r(i)=0;
            c=c+1;
        elseif trm(i)>=t0 %se asigna el valor mayor después de tiempo de corte
            r(i)=max(uz);
        else
            d=d+1; %contador auxiliar
        end
end
for j=1:d %Se construye rampa
        rampa=linspace(0,max(uz),d+1);
        r(c+j)=rampa(j+1);
end
fftrampa=fft(r); %Transformada de rampa
figure, plot(trm,r)
figure, plot(abs(fftrampa))

%Se define la fuente a usar:
F=3;
%=========================================================================

%%
%======================CÁLCULO DE INCOGNITA XI==========================
%Se obtiene incognita para cada frecuencia
for i=1:length(f)
    %Longitud de onda
    lambda=prop.c/f(i);
    dS=lambda/npplo; 
    %Se crean segmentos para cada frecuencia
    [xNf,yNf,nxf,nyf,Sf]=refine_boundarymesh(x,y,min(dS,DS),C);
    
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
    
    for l=1:length(boundary.yCN)
        if boundary.yCN(l)<0
            aux=l+1;
        end
    end
    for k=1:length(boundary.yCN)
        if boundary.yCN(k)==-max(abs(boundary.yCN))
            aux2=k+1;
        end
    end

    %FUENTE
    if F==1 %Desplazamiento horizontal
	    desp=disprig;
    elseif F==2 %Desplazamiento vertical unitario
        disp=[zeros(1,50), ones(1,length(f)-50)];
        uza=zeros(round(length(boundary.xCN)/4-1),length(f));
        uzb=zeros(round(length(boundary.xCN)/4),length(f));
        dispC=[uza;disp;uzb;uzl];
	    desp=dispC(:,i)';
    else   %F==3 Desplazamiento campo cercano
        disp=r(1:length(f));
        uza=zeros(5,length(f));
        uzb=zeros(aux2-5,length(f));
        uzl=zeros(aux-aux2,length(f));
        dispC=[uza;disp;uzb;uzl];
        desp=dispC(:,i)';
    end
   
    %Se construye sistema de ecuaciones y obtiene solución xi
    solution(i).Xi=build_solve_ie(boundary,prop,freq,desp,aux,C);    
    display('freq')
    f(i)
end

%Solución frecuencia "0"
f0=0.0001; 
freq.f=f0;
freq.omega=2*pi*f0;
if F==1
	   desp=disprig;
elseif F==2
	   desp=dispC(:,1)';
else
       desp=dispC(:,1)';
end
solution0.Xi=build_solve_ie(boundary,prop,freq,desp,aux,C);
%========================================================================

%%
%========================================================================
%Definición de puntos de observación
if F==0 %Validación del método con análisis de presiones
    obs=observerl;
    opc=0;
else
    opc=1;
    obs=observer;
end
%========================================================================

%%
%=======================CÁLCULO DESPLAZAMIENTO===========================
for i=1:length(f)
    %for j=1:length(yM) %Para análisis de propagación
    lambda=prop.c/f(i);
    dS=lambda/npplo;    %diferencial para cada longitud de onda
    [xNf,yNf,nxf,nyf,Sf]=refine_boundarymesh(x,y,min(dS,DS),C);  
    
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
    %obs.yO=yy(j,:); %Puntos de observación para malla
    %obs.xO=xx(j,:);
    %Se obtiene desplazamiento y presiones en puntos de observación
    [un(i).p un(i).dispn]=compute_unormal(boundary,prop,freq,solution,obs,C);
    %[unM(i).p(j,:) unM(i).dispn(j,:)]=compute_unormal(boundary,prop,freq,solution,obs,C);
    %end
end

%Calculo de desplazamiento y presión en frecuencia "0"
freq.f=f0;
freq.omega=2*pi*f0;
freq.i=1;
[un0.p un0.dispn]=compute_unormal(boundary,prop,freq,solution0,obs,C);
% for j=1:length(yM) %Frecuencia 0 para malla
%     obs.yO=yy(j,:);
%     obs.xO=xx(j,:);
%     [un0M.p(j,:) un0M.dispn(j,:)]=compute_unormal(boundary,prop,freq,solution0,obs,C);
% end
%========================================================================

%%
%===================TRANSFORMANDO A DOMINIO DEL TIEMPO===================
%CONJUGADO DE RESULTADOS:
    %LÍNEA DE RECEPTORES
%------------------------------PRESIONES    
CC=zeros(length(f)*2,length(obs.xO));
CC(1,:)=real(un0.p(:,1));
for j=1:length(obs.xO)
    for i=2:length(CC(:,1))/2+1
        CC(i,j)=un(i-1).p(j);
        CC(end-i+2,j)=conj(un(i-1).p(j));
    end
end

%--------------------------DESPLAZAMIENTOS    
DD=zeros(length(f)*2,length(obs.xO));
DD(1,:)=real(un0.dispn(:,1));
for j=1:length(obs.xO)
    for i=2:length(DD(:,1))/2+1
        DD(i,j)=un(i-1).dispn(j);
        DD(end-i+2,j)=conj(un(i-1).dispn(j));
    end
end

%TRANSFORMANDO A DOMINIO DEL TIEMPO
for i=1:length(obs.xO)
    dd(:,i)=real(ifft(DD(:,i)));
    cc(:,i)=real(ifft(CC(:,i)));
end

%Definición de vector de tiempo
t=linspace(0,tt,length(f)*2);
dt=t(2)-t(1); %Diferencial de tiempo
for i=1:length(obs.xO)
    ddF(:,i)=computebandfftfilter(dd(:,i),dt,.20,4,.075);
     ddF(:,i)=computebandfftfilter(ddF(:,i),dt,.20,4,.075);
     ccF(:,i)=computebandfftfilter(cc(:,i),dt,.20,4,.075);
     ccF(:,i)=computebandfftfilter(ccF(:,i),dt,.20,4,.075);
   
end
%========================================================================

%%
%==============================FIGURAS===================================
scrsz =get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)]);

for i=1:length(obs.xO)
    plot(t,ddF(:,i)/max(ddF(:,i))+i,'k')
    if M==1
        xlim([0 300])
    else
        xlim([0 900])
    end
    hold on
end
xlabel('Tiempo [s]', 'FontSize',13)
ylabel('Distancia [km]','FontSize',13)
if C==0
    title('Método de imágenes','FontSize',13)
elseif C==1
    title('Condición P=0 en superficie','FontSize',13)
else 
    title('Condición de gravedad','FontSize',13)
end
yticks([0:20:100])
yticklabels({num2str(round(obs.xO(1)/1000)),num2str(round(obs.xO(20)/1000)),num2str(round(obs.xO(40)/1000)),num2str(round(obs.xO(60)/1000)),num2str(round(obs.xO(80)/1000)),num2str(round(obs.xO(100)/1000))})


%Altura máxima:
for i=1:length(obs.xO)
    AM(1,i)=max(dd(:,i));
end
AMAX=max(AM)
%========================================================================
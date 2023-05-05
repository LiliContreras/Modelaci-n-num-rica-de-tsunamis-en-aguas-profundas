%======================================================================
%	  Modelación numérica de tsunamis en aguas profundas
% Aplicación de la teoría de ecuaciones integrales para un modelo de
% aguas profundas.

% Lilibeth Zaira Contreras Alvarado
% Leonardo Ramírez Guzmán
%======================================================================
clear all
close all
clc

%MODELO 1=1     MODELO 2=2
M=1;


%=========================PARÁMETROS DEL MÓDELO=========================
% Definición de frecuencias
frfinal=0.22;%0.22 0.125, 0.15
frinic=0.01; %0.06

if M==1
    tt=600;%300;
else
    tt=800; %Tiempo total de simulación [s] 216
end
w0=2*pi/tt;
df=w0/(2*pi);
f=[frinic:df:frfinal];
%Número de puntos por segmento
npplo=4;

%velocidad 
c=1500; %1500
prop.c=c;
%densidad del agua
rho=1000;
prop.rho=rho;
%número de puntos de observación
nest=100;

if M==1
    %Batimetría Modelo 1:
    x=[ 0 20 160 160.1]*1e3; %metros
    y=[-5 -5 -5  0]*1e3; %metros
else
    %Batimetria Modelo 2:
    x=[0 20 300 700 700.1]*1e3;
    y=[-5 -5 -5 0 0]*1e3;
end

%Se crea malla de batimetría
%	Condiciones:
% Imágenes=0; PresionCero=1; Gravedad=2;
C=0;
DS=1e3;
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

%Grafica de segmentos
figure;
plot(x,y)
hold on
plot(boundary.xN,boundary.yN,'r.','LineStyle','none')
legend('Contorno del modelo','Separación de segmentos','FontSize', 13)
ylabel('Profundidad [m]','FontSize', 14)
xlabel('Distancia x [m]','FontSize', 14)

%variable auxiliar, da índice de donde comienza superficie del modelo
for i=1:length(boundary.yCN)
    if boundary.yCN(i)<0
        aux=i+1;
    end
end
for i=1:length(boundary.yCN)
    if boundary.yCN(i)==-max(abs(boundary.yCN))
        aux2=i+1;
    end
end

%Puntos de observación en superficie
if M==1
    %Modelo 1
    observer.xO=linspace(0,160.1,100)*1e3; %metros
    %observer.yO=0*observer.xO;
    observer.yO=ones(1,length(observer.xO))*(0);
else
    %Modelo 2
    observer.xO=linspace(0,1000,nest)*1e3; %metros
    observer.yO=0*observer.xO;
end

%Puntos de observación en contacto lateral
yO=[0:(yN(aux)-yN(aux-1))/2:5000];
observerl.xO=ones(1,length(yO))*16e4;
observerl.yO=-yO;

%Malla de observación
xM=linspace(0,xN(aux2),aux2);
yM=linspace(0,yN(aux2),10);
[xx yy]=meshgrid(xM,yM);

%=======================================================================


%================================FUENTE=================================
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
rt=0; %50 segundos retraso de tiempo
t0=100; %150 80 Tiempo de corte
trm=linspace(0,tt,length(f)*2);
r=ones(1,length(f)*2);
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
for j=1:d
        rampa=linspace(0,max(uz),d+1);
        r(c+j)=rampa(j+1);
end
fftrampa=fft(r);
figure, plot(trm,r)
figure, plot(abs(fftrampa))
uza=zeros(5,length(f));
uzb=zeros(aux2-5,length(f));
uzl=zeros(aux-aux2,length(f));
R=[uza;fftrampa(1:length(f));uzb;uzl];

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
    %disp=ones(1,length(f));
    %%disp=fftrampa(1:length(f));
    %uza=zeros(round(length(boundary.xCN)/4-1),length(f));
    %uzb=zeros(round(length(boundary.xCN)/4),length(f));
    %disp=[uza;disp;uzb;uzl];
    %%uza=zeros(5,length(f));
    %%uzb=zeros(aux2-5,length(f));
    %%uzl=zeros(aux-aux2,length(f));
    %%disp=[uza;disp;uzb;uzl];

    %FUENTE
    if F==1
	    desp=disprig;
    elseif F==2
        disp=ones(1,length(f));
        uza=zeros(round(length(boundary.xCN)/4-1),length(f));
        uzb=zeros(round(length(boundary.xCN)/4),length(f));
        disp=[uza;disp;uzb;uzl];
	    desp=disp(:,i)';
    else
	    %desp=Fuz(i,:);
        disp=fftrampa(1:length(f));
        uza=zeros(5,length(f));
        uzb=zeros(aux2-5,length(f));
        uzl=zeros(aux-aux2,length(f));
        disp=[uza;disp;uzb;uzl];
        desp=disp(:,i)';
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
	   desp=disp(:,1)';
else
	   %desp=Fuz(1,:);
       desp=R(:,1)';
end
solution0.Xi=build_solve_ie(boundary,prop,freq,desp,aux,C);
%========================================================================
%{
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
%}
opc=0;
obs=observer;

%Ricker:
fm=0.05; %0.2,0.1
t0=rt; 
t1=linspace(0,tt,length(f)*2);
w=(1-2*(pi^2)*(fm^2).*(t1-t0).^2).*exp(-(pi^2)*(fm^2).*(t1-t0).^2);

W=fft(w);
figure, plot(t1,w)
figure, plot(abs(W))

%Filtro:
for i=1:length(f)*2
    WP(i)=W(i);
    if W(i)<0.6
        WP(i)=0;
    end
end



%=======================CÁLCULO DESPLAZAMIENTO===========================
for i=1:length(f)
    %for j=1:length(yM)
    lambda=prop.c/f(i);
    dS=lambda/npplo;    
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
    %obs.yO=yy(j,:);
    %obs.xO=xx(j,:);
    %Se obtiene desplazamiento y presiones en puntos de observación
    [un(i).p un(i).dispn]=compute_unormal(boundary,prop,freq,solution,obs,C);
    %[unM(i).p(j,:) unM(i).dispn(j,:)]=compute_unormal(boundary,prop,freq,solution,obs,C);
    %end
end

%Calculo de frecuencia "0"
freq.f=f0;
freq.omega=2*pi*f0;
freq.i=1;
[un0.p un0.dispn]=compute_unormal(boundary,prop,freq,solution0,obs,C);
%for j=1:length(yM)
%    obs.yO=yy(j,:);
%    obs.xO=xx(j,:);
%    [un0M.p(j,:) un0M.dispn(j,:)]=compute_unormal(boundary,prop,freq,solution0,obs,C);
%end
%========================================================================

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
%CONVOLUCIÓN CON ONDÍCULA
for i=1:length(obs.xO)
    for j=1:length(f)*2
        CCR(j,i)=CC(j,i).*W(j);
    end
end
%TRANSFORMANDO A DOMINIO DEL TIEMPO
for i=1:length(obs.xO)
    cc(:,i)=real(ifft(CCR(:,i)));
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
%CONVOLUCIÓN CON ONDÍCULA
for i=1:length(obs.xO)
    for j=1:length(f)*2
        DDR(j,i)=DD(j,i)*W(j);
    end
end
%TRANSFORMANDO A DOMINIO DEL TIEMPO
for i=1:length(obs.xO)
    dd(:,i)=real(ifft(DDR(:,i)));
end
%=======================================================================

%{
%MALLA DE RECEPTORES:
%------------------------------PRESIONES  
%CCM=(x,frecuencias)
CCM=


CCM=zeros(length(yM),length(xM),length(f)*2);
CCM(:,:,1)=real(un0M.p(:,:));
for k=1:length(yM)
    for j=1:length(xM)
        for i=2:length(CCM(1,1,:))/2+1
            CCM(k,j,i)=unM(i-1).p(k,j);
            CCM(k,j,end-i+2)=conj(unM(i-1).p(k,j));
        end
    end
end
%CONVOLUCIÓN CON ONDÍCULA
for j=1:length(yM)
    for i=1:length(xM)
        for k=1:length(f)*2
            CCMW(j,i,k)=CCM(j,i,k)*W(k);
        end
    end
end
%TRANSFORMANDO A DOMINIO DEL TIEMPO
for j=1:length(yM)
    for i=1:length(xM)
        ccm(j,i,:)=real(ifft(CCMW(j,i,:)));
    end    
end



for i=1:length(f)
   for j=1:length(obs.xO)
      Ad(i,j)=un(i).dispn(j); 
      Ap(i,j)=un(i).p(j); 
   end
end

for i=1:length(f)
   for j=1:length(yM)
      AdM(j,:,i)=unM(i).dispn(j,:); 
      ApM(j,:,i)=unM(i).p(j,:); 
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
%}
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
    tmax=max(x)/prop.c;
    figure
    view(3)
    for i=1:length(obs.xO)
        xo=ones(1,length(f)*2)*obs.xO(i);
        t=linspace(0,tt,length(f)*2);
        hold on
        plot3(xo,t,dd(:,i),'b')
        zlim([-2 2])
    end
    title('Desplazamientos normales')
    xlabel('Estaciones')
    ylabel('Tiempo [s]')
    zlabel('Altura de ola')
    ylim([0 300])
    figure
    view(3)
    for i=1:length(obs.xO)
        xo=ones(1,length(f)*2)*obs.xO(i);
        t=linspace(0,tt,length(f)*2);
        hold on
        plot3(xo,t,cc(:,i),'b')
        %zlim([-3 3])
    end
    title('Presiones-Imágenes')
    xlabel('Estaciones')
    ylabel('Tiempo [s]')
    zlabel('Presión')
end
%{
%Mallas
figure
view(2)
for i=1:length(f)*2
    surf(xx,yy,ccm(:,:,i))
    shading interp
    clim([-100 100])
    view(2)
    pause(1)
end
figure
%view(2)
for i=1:length(f)
    surf(xx,yy,UNM(:,:,i))
    shading interp
    %view(2)
    pause(0.3)
end
%}
%for i=1:length(xM)
%    CCM1(i,:)=ccm(1,i,:);
%end

figure
for i=1:length(obs.xO)
    subplot(3,1,1)
    plot(abs(W))
    subplot(3,1,2)
    plot(abs(DD(:,i)))
    ylim([0 60])
    %pause(0.3)
    subplot(3,1,3)
    plot(abs(DDR(:,i)))
    ylim([0 60])
    pause(0.3)
end

for i=1:length(obs.xO)
    DDT(:,i)=DD(:,i)*tukeywin(length(f)*2,0.1);
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

%Ciclo de estaciones
t=linspace(0,tt,length(f)*2);
scrsz =get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)]);
for i=1:length(obs.xO)
    plot(t,dd(:,i)/max(dd(:,i))+i,'k')
    hold on
end
xlabel('Tiempo [s]', 'FontSize',13)
ylabel('Estaciones','FontSize',13)
if C==0
    title('Método de imágenes','FontSize',13)
elseif C==1
    title('Condición P=0 en superficie','FontSize',13)
else 
    title('Condición de gravedad','FontSize',13)
end
%xlim([0 300])

%ANIMACIÓN PROPAGACIÓN MALLA
 figure
 filename = 'Propagacion_Tiempo600s.gif';
view(2)
numberOfFrames=length(f)*2;
%[x, y] = meshgrid(x1d, y1d);
for frameIndex = 1 : numberOfFrames
	z = ccm(:,:,frameIndex);
	%cla reset;
	surf(xx,yy,z);
    shading interp
    view(2)
	%axis('tight')
	%zlim([-1000, 1000]);
	clim([-150 150])
    %axis equal
    %xlim([73319 111486]);
    axis equal
    drawnow;
    pause(1)
	%thisFrame = getframe(gca);
	% Write this frame out to a new video file.
%  	writeVideo(writerObj, thisFrame);
	%myMovie(frameIndex) = thisFrame;
    frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if frameIndex == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end
end


%FIGURA DE FUENTE:
[yf,xf]=meshgrid(f,xN(1:aux+1));

figure
subplot(2,2,1)
plot(xN(1:aux),uz)
ylabel('Desplazamiento','FontSize',13)
subplot(2,2,2)
plot(x,y,LineWidth=2)
hold on
plot(boundary.xN,boundary.yN,'r.','LineStyle','none')
ylabel('Profundidad [m]','FontSize',13)
subplot(2,2,3)
plot(t,r)
ylabel('Desplazamiento','FontSize',13)
xlabel('Tiempo [s]','FontSize',13)
subplot(2,2,4)
view(2)
surf(xf,yf,abs(R)); shading interp
ylabel('Frecuencia [Hz]','FontSize',13)
xlabel('Distancia [m]','FontSize',13)
%zlim([0 200])
colorbar



%
%  Se construye y resuelve el sistema de ecuaciones integrales
%

function Xi=build_solve_ie(boundary,prop,freq,disprig,aux,condicion)
%Se construye el sistema según la condición
if condicion==0 %Método de imágenes
					% Matriz G
   for iC=1:length(boundary.xCN)/2            
       piC.x=boundary.xCN(iC); %Puntos centrales de segmentos  
       piC.y=boundary.yCN(iC);
       n.n  =[boundary.nxN(iC) boundary.nyN(iC)];
       for iS=1:length(boundary.xCN)/2
       	   %Segmentos de modelo original
           iA=iS;
           iB=iA+1;
           sourceregP.x=[boundary.xN(iA) boundary.xN(iB)];
           sourceregP.y=[boundary.yN(iA) boundary.yN(iB)];
           %Segmentos de imágen
           iA=length(boundary.xN)-iS;
           iB=iA+1;
           sourceregN.x=[boundary.xN(iA) boundary.xN(iB)];
           sourceregN.y=[boundary.yN(iA) boundary.yN(iB)];

       	   if(iC==iS)     %Cálculo de contribución sobre el mismo segmento      
             A(iC,iS)=   -.5    -inte_ic_is(piC,n,sourceregN,prop,freq);
           else  %Calculo de contribución de otros segmentos
             A(iC,iS)=inte_ic_is(piC,n,sourceregP,prop,freq)-inte_ic_is(piC,n,sourceregN,prop,freq);
           end        
      end
  end
					% Matriz U
  for iC=1:length(boundary.xCN)/2   
      n  =[boundary.nxN(iC) boundary.nyN(iC)];
      if length(disprig)>2
         nuv=[0 disprig(iC)];
         B(iC,1)=dot(n,nuv);
      else	
         B(iC,1)=dot(n,disprig);   % Desplazamiento vertical
      end
  end

  % Se invierte Xi
  Xi=A\B;

elseif condicion==1 %Condición P=0 en superficie
% 				Matriz G
%  Fondo
       for iC=1:aux       %segmentos del fondo del modelo    
   	   piC.x=boundary.xCN(iC);  
   	   piC.y=boundary.yCN(iC);
   	   n.n  =[boundary.nxN(iC) boundary.nyN(iC)];
   	   for iS=1:length(boundary.xCN) %POSICIÓN SEGMENTO
       		%Segmento de contribución en el fondo
       	       iA=iS; 
       	       iB=iA+1;
       	       sourceregP.x=[boundary.xN(iA) boundary.xN(iB)];
               sourceregP.y=[boundary.yN(iA) boundary.yN(iB)];
       	       if(iC==iS)   %Contribución sobre el mismo segmento                         
                  A(iC,iS)= -.5;
               else  %Contribución de los otros segmentos
                  A(iC,iS)=inte_ic_is(piC,n,sourceregP,prop,freq);
               end      
            end
       end
%  Superficie
       for iC=aux:length(boundary.xCN)     %segmentos en superficie      
   	   piC.x=boundary.xCN(iC);  
   	   piC.y=boundary.yCN(iC);
   	   n.n  =[boundary.nxN(iC) boundary.nyN(iC)];
   	   for iS=1:length(boundary.xCN) %POSICIÓN SEGMENTO
       
       	       iA=iS; 
       	       iB=iA+1;
       	       sourceregP.x=[boundary.xN(iA) boundary.xN(iB)];
       	       sourceregP.y=[boundary.yN(iA) boundary.yN(iB)];
               if(iC==iS)                               
                  A(iC,iS)= inteP_ic_is(piC,sourceregP,prop,freq);%-.5;
               else
                  A(iC,iS)=inteP_ic_is(piC,sourceregP,prop,freq);
               end      
           end
       end

% 				Matrix U

       for iC=1:aux   %Desplazamiento en segmentos del fondo
    	   n  =[boundary.nxN(iC) boundary.nyN(iC)];
	   if length(disprig)>2
	      nuv=[0 disprig(iC)];
	      B(iC,1)=dot(n,nuv);
       else
    	      B(iC,1)=dot(n,disprig);   % Desplazamiento vertical
	   end
       end
       for iC=aux+1:length(boundary.xCN) %Desplazamiento en segmentos de superficie
           B(iC,1)=0;
       end
% Se invierte Xi
Xi=A\B;

else %Condición de gravedad
% 				Matriz G
%  Fondo
    for iC=1:aux           
        piC.x=boundary.xCN(iC);  
        piC.y=boundary.yCN(iC);
        n.n  =[boundary.nxN(iC) boundary.nyN(iC)];
        for iS=1:length(boundary.xCN) %POSICIÓN SEGMENTO
       
            iA=iS; 
            iB=iA+1;
            sourceregP.x=[boundary.xN(iA) boundary.xN(iB)];
            sourceregP.y=[boundary.yN(iA) boundary.yN(iB)];
       
            if(iC==iS)                               
               A(iC,iS)= -.5;
            else
               A(iC,iS)=inte_ic_is(piC,n,sourceregP,prop,freq);
            end      
    
        end
    end
%  Superficie
    g=9.81;
    for iC=aux:length(boundary.xCN)           
        piC.x=boundary.xCN(iC);  
        piC.y=boundary.yCN(iC);
        n.n  =[boundary.nxN(iC) boundary.nyN(iC)];
        for iS=1:length(boundary.xCN) %POSICIÓN SEGMENTO
       
            iA=iS; 
            iB=iA+1;
            sourceregP.x=[boundary.xN(iA) boundary.xN(iB)];
            sourceregP.y=[boundary.yN(iA) boundary.yN(iB)];

            if(iC==iS)                               
               A(iC,iS)= (freq.omega/g)*inteP_ic_is(piC,sourceregP,prop,freq)+inte_ic_is(piC,n,sourceregP,prop,freq);
            else
            A(iC,iS)=(freq.omega/g)*inteP_ic_is(piC,sourceregP,prop,freq)+inte_ic_is(piC,n,sourceregP,prop,freq);
            end      
    
        end
    end

% 				Matriz U

    for iC=1:aux   
        n  =[boundary.nxN(iC) boundary.nyN(iC)];
	if length(disprig)>2
           nuv=[0 disprig(1,iC)];
	   B(iC,1)=dot(n,nuv);
	else
           B(iC,1)=dot(n,disprig);   % Desplazamiento vertical
	end
    end
    for iC=aux+1:length(boundary.xCN)
        B(iC,1)=0;
    end
% Se invierte Xi
Xi=A\B;
end

end

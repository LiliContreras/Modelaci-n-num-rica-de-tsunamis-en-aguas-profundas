%
%  build and solve the equation
%

function [p ,dispn]=compute_unormal(boundary,prop,freq,solution,observer,condicion)    

if condicion==0 %Método de imágenes
    %                    Matriz G
    for iO=1:length(observer.xO)            
        piO.x=observer.xO(iO);  
        piO.y=observer.yO(iO);
        n.n  =[0 1];
        for iS=1:length(boundary.xCN)/2
       
            iA=iS;
            iB=iA+1;
            sourceregP.x=[boundary.xN(iA) boundary.xN(iB)];
            sourceregP.y=[boundary.yN(iA) boundary.yN(iB)];
       
            iA=length(boundary.xN)-iS;
            iB=iA+1;
            sourceregN.x=[boundary.xN(iA) boundary.xN(iB)];
            sourceregN.y=[boundary.yN(iA) boundary.yN(iB)];
       
            AP(iO,iS)=inteP_ic_is(piO,sourceregP,prop,freq)-inteP_ic_is(piO,sourceregN,prop,freq);
            AUN(iO,iS)=inte_ic_is(piO,n,sourceregP,prop,freq)-inte_ic_is(piO,n,sourceregN,prop,freq);
       
        end
    end
    p=AP*solution(freq.i).Xi;
    dispn=AUN*solution(freq.i).Xi;
else %Condición P=0 en superficie o gravedad
    for iO=1:length(observer.xO)            
        piO.x=observer.xO(iO);  
        piO.y=observer.yO(iO);
        n.n=[0 1];
        for iS=1:length(boundary.xCN)
       
            iA=iS;
            iB=iA+1;
            sourceregP.x=[boundary.xN(iA) boundary.xN(iB)];
            sourceregP.y=[boundary.yN(iA) boundary.yN(iB)];
       
            AP(iO,iS)=inteP_ic_is(piO,sourceregP,prop,freq);
            AUN(iO,iS)=inte_ic_is(piO,n,sourceregP,prop,freq);
       
        end
    end
end
p=AP*solution(freq.i).Xi;
dispn=AUN*solution(freq.i).Xi;

end

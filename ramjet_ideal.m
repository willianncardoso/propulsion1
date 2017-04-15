%zerando------------------------------------------------------------------%
clear variables
close all
clc
format long
%-------------------------------------------------------------------------%

% Dados iniciais conhecidos-----------------------------------------------%
T0=216.7;                   %temperatura inicial    |K|
y = 1.4;                    %                       |J/kgK|
cp = 1004;                  %                       |kJ/kg| 
hpr = 42800;                              
Tt4 = [1600; 1900; 2200];   %temperaturas no fim da camara de combustao |K| CONJUNTO FECHADO()
M0 = 0:0.1:7;

j = length(Tt4);%cada linha do gr?fico CONJUNTO FECHADO
n = length(M0);
%-------------------------------------------------------------------------%

%saida--------------------------------------------------------------------%
%       empuxo = |N/(kg/s)|
%       f      = |adimensional| razao combustivel/ air
%       s      = |(mg/s)/N| consumo espeficico de combustivel por empuxo
%       nT     = |adimensional| eficiencia termica
%       nP     = |adimensional| eficiencia propulsiva
%       n0     = |adimensional| eficiencia total
%-------------------------------------------------------------------------%

%pre-alocando matrizes com zeros para preenche-las------------------------%
tal_r   = zeros(j,n); 
tal_lambida = zeros(j,n); 
raz_1   = zeros(j,n);
raz_2   = zeros(j,n);
f       = zeros(j,n);
s       = zeros(j,n);
nP      = zeros(j,n);
nT      = zeros(j,n);
n0      = zeros(j,n);
%-------------------------------------------------------------------------%

%loop de varredura em M0--------------------------------------------------%
for i = 1:j   
    for k = 1:n
        
        R=((y-1)/y)*cp;
        a0=sqrt(y*R*T0);
        tal_r=1+((y-1)/2).*(M0.^2);
        tal_lambida=(Tt4)/(T0);
        raz_1=M0.*sqrt(tal_lambida./tal_r);
        %raz_1 = V9/a0
        raz_2=a0.*((raz_1)-M0);
        %raz_2 = F/m0ponto = Empuxo motor n?o instalado
        f=((cp*T0)/hpr).*(tal_lambida-tal_r)/1000;
        s=(f./(raz_2)).*1000000;
        nT=1-(1./tal_r);
        nP=2./(sqrt(tal_lambida./tal_r)+1);
        n0=(2.*(tal_r-1))./(sqrt(tal_lambida.*tal_r)+tal_r);
        
        %Graficos---------------------------------------------------------%
        figure(1)
        subplot(2,2,1)
        plot(M0,raz_2)
        xlabel('M0')
        ylabel('empuxo')
        
        subplot(2,2,2)
        plot(M0,f)
        xlabel('M0')
        ylabel('f')
        
        subplot(2,2,3)
        plot(M0,s)
        xlabel('M0')
        ylabel('s')
        
        subplot(2,2,4)
        plot(M0,nT,'--',M0,nP,M0,n0)
       
        xlabel('M0')
        ylabel('n')
    end
    
end  
%-------------------------------------------------------------------------%

%Mostrando em linha de comando--------------------------------------------%

%-------------------------------------------------------------------------%




   


% Turbofan ideal 
% aqui corremos em pi_c

%zerando------------------------------------------------------------------%
clear variables
close all
clc
format long
%-------------------------------------------------------------------------%

% Dados iniciais conhecidos-----------------------------------------------%
M0   = 0.9;                     %M0
T0   = 216.7;                   %temperatura inicial    |K|
y    = 1.4;                     %                       |J/kgK|
cp   = 1004;                    %                       |kJ/kg| 
hpr  = 42800;                              
Tt4  = 1670;                    %temperaturas no fim da camara de combustao |K|
pi_c = 10:1:40;                 %varredura em pi_c
PI_f = 2:1:5;                   %CONJUNTO FECHADO
pi_f = PI_f';                   %transposta 
                

j = length(pi_f);%cada linha do grafico CONJUNTO FECHADO
n = length(pi_c);%onde corremos   
%-------------------------------------------------------------------------%

%saida--------------------------------------------------------------------%
%       empuxo = |N/(kg/s)|
%       f      = |adimensional| razao combustivel/ air
%       f0     = 
%       s      = |(mg/s)/N| consumo espef?cico de combust?vel por empuxo
%       nT     = |adimensional| eficiencia termica
%       nP     = |adimensional| eficiencia propulsiva
%       n0     = |adimensional| eficiencia total
%       pi_f ou alfa
%-------------------------------------------------------------------------%

%pre-alocando matrizes com zeros para preenche-las------------------------%
tal_r   = zeros(j,n); 
tal_lambida = zeros(j,n); 
raz_1   = zeros(j,n);
raz_2   = zeros(j,n);
ftotal  = zeros(j,n);
s       = zeros(j,n);
nP      = zeros(j,n);
nT      = zeros(j,n);
n0      = zeros(j,n);
%-------------------------------------------------------------------------%



%loop de varredura em pi_c--------------------------------------------------%
for i = 1:j    
    for k = 1:n

    R           = ((y-1)/y)*cp;
    a0          = sqrt(y*R*T0);
    tal_r       = 1+((y-1)/2).*(M0.^2);
    tal_lambida = (Tt4)/(T0);
    tal_c       = (pi_c).^((y-1)/y);
    
        if pi_f     ~= 0 %pi_f especificado
            tal_f       = (pi_f).^((y-1)/y);  
            alfa        = ((tal_lambida.*(tal_c-tal_f))./((tal_r.*tal_c).*(tal_f-1)))-((tal_c-1)./(tal_f-1));
            if alfa < 0
                alfa = 0;
            else
            end
        else             %alfa especificado
            tal_f       = ((tal_lambida./tal_r)-(tal_c-1)+alfa)./((tal_lambida/(tal_r.*tal_c))+alfa);
            pi_f        = (tal_f)^(y/(y-1));
        end
        
    tal_t       = 1-(tal_r./tal_lambida).*(tal_c-1+alfa.*(tal_f-1));
    f           = ((cp*T0)/hpr).*(tal_lambida-tal_r.*tal_c)/1000;
    tal_M       = (1./(1+alfa)).*(1+alfa.*((tal_r.*tal_f)./(tal_lambida.*tal_t)));
    %raz_4 = T9/T0
    raz_4       = (tal_lambida.*tal_t.*tal_M)./(tal_r.*tal_f);
    M9          = sqrt(((2/(y-1)).*(tal_r.*tal_f-1)));
    %raz_1 = V9/a0
    raz_1       = sqrt(raz_4).*M9;
    f0          = f./(1+alfa);
    %raz_2 = F/m0ponto EMPUXO
    raz_2       = a0.*(raz_1-M0);
    s           = (f0./raz_2).*1000000;
    nT          = (((y-1)/2).*((cp.*T0)./(f0.*hpr)).*(raz_1.^2-M0.^2));
    nP          = (2.*M0)./(raz_1+M0);
    n0          = (nT.*nP);
  
    
    %Graficos---------------------------------------------------------%
        figure(1)
        subplot(2,2,1)
        plot(pi_c,raz_2)
        xlabel('pi_c')
        ylabel('empuxo')
        grid
        legend('pi_f1','pi_f2','pi_f3','pi_f4')
        
        subplot(2,2,2)
        plot(pi_c,s)
        xlabel('pi_c')
        ylabel('s')
        grid
       
        subplot(2,2,3)
        plot(pi_c,alfa)
        xlabel('pi_c')
        ylabel('alfa')
        grid
       

    end
end

        
     

%----------------------------------------------------------------------%


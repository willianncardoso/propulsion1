% Turbofan ideal 
% aqui corremos em pi_c

%zerando------------------------------------------------------------------%
clear variables
close all
clc
format long
%-------------------------------------------------------------------------%

% Dados iniciais conhecidos-----------------------------------------------%
M0   = 0.8;                     %M0
T0   = 240;                   %temperatura inicial    |K|
y    = 1.4;                     
cp   = 1004;                    %                       |J/kgK| 
hpr  = 42800;                   %                       |kJ/kg|           
Tt4  = 1370;                    %temperaturas no fim da camara de combustao |K|
nprop = 0.83;
pi_c = 0:1:40;                 %varredura em pi_c
tal_T = 0.4:0.1:0.8;           %CONJUNTO FECHADO
tal_t = tal_T';           
                

j = length(tal_t);%cada linha do gr?fico CONJUNTO FECHADO
n = length(pi_c);%onde corremos   
%-------------------------------------------------------------------------%

%saida--------------------------------------------------------------------%
%       empuxo = |N/(kg/s)|
%       f      = |adimensional| razao combustivel/ air
%       f0     = 
%       s      = |(mg/s)/N| consumo espef?cico de combust?vel por empuxo
%       nT     = |adimensional|
%       nP     = |adimensional|
%       n0     = |adimensional|
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

R           = ((y-1)/y)*cp;
    a0          = sqrt(y*R*T0);
    tal_r       = 1+((y-1)/2).*(M0.^2);
    tal_lambida = (Tt4)/(T0);


%loop de varredura em pi_c--------------------------------------------------%
for i = 1:j    
    for k = 1:n

   
    tal_c       = (pi_c).^((y-1)/y);
    f           = ((cp*T0)/hpr).*(tal_lambida-tal_r.*tal_c)/1000;
    tal_tH      = 1 - (tal_r/tal_lambida).*(tal_c-1);
    tal_tL      = tal_t./tal_tH;
    %raz_1 = V9/a0
    raz_1       = sqrt((2/(y-1)).*((tal_lambida.*tal_t)-((tal_lambida)./(tal_r.*tal_c))));
    Cc          = (y-1).*M0.*(raz_1-M0);
    Cprop       = nprop.*tal_lambida.*tal_tH.*(1-tal_tL);
    Ctot        = Cprop+Cc;
    %raz_2 = F/m0ponto EMPUXO
    raz_2       = (Ctot.*cp.*T0)./(a0.*M0);
    s           = (f./raz_2).*1000000;
    nT          = (1-(1./(tal_r.*tal_c)));
    n0          = (Ctot)./(tal_lambida-tal_r.*tal_c);
    nP          = (n0./nT);
  
    
    %Graficos---------------------------------------------------------%
        figure(1)
        subplot(2,2,1)
        plot(pi_c,raz_2)
        xlabel('pi_c')
        ylabel('empuxo')
        grid
        legend('1','2','3','4')
        
        subplot(2,2,2)
        plot(pi_c,s)
        xlabel('pi_c')
        ylabel('s')
        grid
       
        subplot(2,2,3)
        plot(pi_c,tal_t)
        xlabel('pi_c')
        ylabel('tal_t')
        grid
        
        disp(raz_2)

    end
end

        
     

%----------------------------------------------------------------------%


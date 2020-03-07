%Lumped parameter thermal network

%% Initilization
% Parameter table
% Ns       - stator poles
% Nr       - rotor poles
% Arc_slot - slot arc 
% Arc_st   - Stator tooth arc
% Ro       - Stator slot outter radill mm
% Ri       - Stator slot inner radill
% Rs       - Stator outter surface
% g        - Air gap length
% Arc_rt   - Rotor tooth arc
% ro       - Rotor outter radill
% ri       - Rotor inner radill
% r_shaft  - Rotor shaft radill
% Q        - slot number
% hc       - copper thermal conductivity
% hl       - liner thermal conductivity
% hr       - resin thermal conductivity
% hs       - steel thermal conductivity
% ho       - outer air convection conductivity
% hi       - inner air convection conductivity
% hh       - housing thermal conductivity
% hair     - air thermal conductivity
% kc       - copper packing factor
% kl       - liner packing factor
% kr       - resin packing factor
% cc       - copper heat capacity
% csteel   - stator steel heat capacity
% cair     - air heat capacity
% mc       - copper mass
% mst      - stator tooth mass
% mstb     - stator tooth back iron mass
% mssb     - stator slot back iron mass
% mairgap  - airgap mass
% mairrs   - rotor slot air mass
% mrt      - rotor tooth mass
% mrb      - rotor back iron mass
% Rh       - Housing radill
% ch       - Housing heat capacity
% mh       - Housing mass
% delta    - Time step


Ns = 12;
Nr = 10;
Arc_slot = 7.5*2*pi/180;
Arc_st = 7.5*2*pi/180;
Ro = 90/1000; 
Ri = 60/1000;
Rs = 105/1000;
r_shaft = 25/1000;
g = 0.5/1000;
Arc_rt = 12*2*pi/180;
Arc_rslot = 24*2*pi/180;
ro = 59.5/1000;
ri = 45/1000;
Q  = 24;
Lstk = 120/1000;
hc = 0.2;
hs = 42.5;
hair = 0.0233;
ho = 20;  %GJ Li paper's data 30W/mm^2
hi = 40;
hh =  42.5 ; %refer to hs
kc = 0.5;
Nu = 2;
e = (ro+g/2)*2*pi;
h_cov_air = Nu*2*e/hair;
cc = 390;
csteel = 480;
cair = 1000;
ch  = 480; % refer to steel 
mc = 8.1;
mst = 7;
mstb = 1.582;
mssb = 1.468;
mairgap = 0.022;
mairrs = 0.396;
mrt = 2.06;
mrb = 3.48;
mh = 33.38;
Rh = 140/1000;
delta = 0.18;
FlagCon = 0;

%% Conduction resistance
% Winding
Rcx  = Arc_slot*(Ro+Ri)/2/(hc*kc*(Ro-Ri)*Lstk)/Q;
Rcy  = Arc_slot*log(Ro/Ri)/(hc*kc*Lstk)/Q;
Rix  = 0;
Riy  = 0;
Rrx  = 0;
Rry  = 0;
Rcirx = Rcx + Rix + Rrx;
Rciry = Rcy + Riy + Rry;

% Stator back iron
Rsby_slot = Arc_slot*log(Rs/Ro)/(hs*Lstk)/Q;
Rsby_tooth = Arc_st*log(Rs/Ro)/(hs*Lstk)/Q;
Rsbx = (Arc_slot+Arc_st)/2*(Rs+Ro)/2/(hs*(Rs-Ro)*Lstk)/(2*Q);

% Housing
Rh_tooth = 2*pi*(Arc_st)/Q*log(Rh/Rs)/(hh*Lstk);
Rh_slot  = 2*pi*(Arc_slot)/Q*log(Rh/Rs)/(hh*Lstk);
Rhx = (Arc_slot+Arc_st)/2*(Rs+Rh)/2/(hh*(Rh-Rs)*Lstk)/(2*Q);

% Stator tooth
Rstx = (Ro-Ri)/(hs*Lstk*Arc_st*Ri)/Q/2;
Rsty = (Arc_st*Ri)/(hs*Lstk*(Ro-Ri))/Q;

% Air gap
Rcov_airgap = 1/(h_cov_air*e*Lstk);
Rcon_airgap_rt = Arc_rt*log((ro+g)/ro)/(hair*Lstk)/Nr;
Rcon_airgap_rsloty = Arc_rslot*log((ro+g)/ri)/(hair*Lstk)/Nr;
Rcon_airgap_rslotx = Arc_rslot*(ro+g+ri)/2/(hair*(g+ro-ri)*Lstk)/Nr;

% Rotor tooth
Rcon_rty = Arc_rt*log((ro)/ri)/(hs*Lstk)/Nr;
Rcon_rtx = Arc_rt*(ro+ri)/2/(hs*(ro-ri)*Lstk)/Nr;

% Rotor back iron
Rrby = 2*pi*log(ri/r_shaft)/(hs*Lstk);

% Outter and inner air convection
Rcov_o = 1/(ho*Rh*2*pi*Lstk);
Rcov_i = 1/(hi*ri*2*pi*Lstk);

Rcovo_slot  = (Arc_slot+Arc_st)/(ho*Rh*2*pi*Lstk)/(Arc_slot); 
Rcovo_tooth = (Arc_slot+Arc_st)/(ho*Rh*2*pi*Lstk)/(Arc_st);

% Heat radiation




%% Heat capacity
Cst = csteel*mst;
Ccopper = cc*mc;
Cstb = csteel*mstb;
Cssb = csteel*mssb;
Cairgap = cair*mairgap;
Cairrs = cair*mairrs;
Crt = csteel*mrt;
Crb = csteel*mrb;
Ch_slot = ch*mh*Arc_slot/(Arc_st+Arc_slot);
Ch_tooth = ch*mh*Arc_st/(Arc_st+Arc_slot);
C = [Ccopper Cst Cssb Cstb Cairgap Cairrs Crt Crb Ch_slot Ch_tooth];

%% Power Loss injection

P = [100 0 0 0 0 0 0 0 0 0];



%% Calculation
% EXample case for thermal circuit cal;

%delta = 1e-3; 
%dT = [100 100];
%j = 1;
%T = [0 0];
%while(j<100000)
%Y = [1001 -1;-1 1];
%if j >1
%Q = [1;1]*T(j-1,:).*Y;
%end
%if j == 1
%    Q = [0 0 ;0 0];
%end
%dQ = [Q(1,1)+Q(1,2);Q(2,2)+ Q(2,1)];
%dP = [1000*delta-dQ(1)*delta;-dQ(2)*delta];
%dT = [dP(1)*delta/1 dP(1)*delta/1];
%if j > 1
%T(j,:) = T(j-1,:)+dT;
%else if j == 1
%        T(j,:) = dT+T;
%    end
%end
%j = j+1
%end

% Martrix definition

R(1,1) = ;
R(1,2) = ;
R(1,3) = ;
R(1,4) = ;
R(1,5) = ;
R(1,6) = ;
R(1,7) = ;
R(1,8) = ;
R(1,9) = ;

R(2,1) = ;
R(2,2) = ;
R(2,3) = ;
R(2,4) = ;
R(2,5) = ;
R(2,6) = ;
R(2,7) = ;
R(2,8) = ;
R(2,9) = ;

R(3,1) = ;
R(3,2) = ;
R(3,3) = ;
R(3,4) = ;
R(3,5) = ;
R(3,6) = ;
R(3,7) = ;
R(3,8) = ;
R(3,9) = ;

R(4,1) = ;
R(1,2) = ;
R(1,3) = ;
R(1,4) = ;
R(1,5) = ;
R(1,6) = ;
R(1,7) = ;
R(1,8) = ;
R(1,9) = ;

R(1,1) = ;
R(1,2) = ;
R(1,3) = ;
R(1,4) = ;
R(1,5) = ;
R(1,6) = ;
R(1,7) = ;
R(1,8) = ;
R(1,9) = ;

R(1,1) = ;
R(1,2) = ;
R(1,3) = ;
R(1,4) = ;
R(1,5) = ;
R(1,6) = ;
R(1,7) = ;
R(1,8) = ;
R(1,9) = ;

R(1,1) = ;
R(1,2) = ;
R(1,3) = ;
R(1,4) = ;
R(1,5) = ;
R(1,6) = ;
R(1,7) = ;
R(1,8) = ;
R(1,9) = ;

R(1,1) = ;
R(1,2) = ;
R(1,3) = ;
R(1,4) = ;
R(1,5) = ;
R(1,6) = ;
R(1,7) = ;
R(1,8) = ;
R(1,9) = ;

R(1,1) = ;
R(1,2) = ;
R(1,3) = ;
R(1,4) = ;
R(1,5) = ;
R(1,6) = ;
R(1,7) = ;
R(1,8) = ;
R(1,9) = ;

R(1,1) = ;
R(1,2) = ;
R(1,3) = ;
R(1,4) = ;
R(1,5) = ;
R(1,6) = ;
R(1,7) = ;
R(1,8) = ;
R(1,9) = ;

R(1,1) = ;
R(1,2) = ;
R(1,3) = ;
R(1,4) = ;
R(1,5) = ;
R(1,6) = ;
R(1,7) = ;
R(1,8) = ;
R(1,9) = ;

R(1,1) = ;
R(1,2) = ;
R(1,3) = ;
R(1,4) = ;
R(1,5) = ;
R(1,6) = ;
R(1,7) = ;
R(1,8) = ;
R(1,9) = ;
% Resistance2Conductance
[m m] = size(R);
Y=R;
for index1 =1:1:m
    for index2 =1:1:m
        if(R(index1,index2)~=0)
            Y(index1,index2) = 1/R(index1,index2);
        else if(1)
                Y(index1,index2)=0;
            end
        end
    end
end
for index1 = 1:1:m
    Y(index1,index1) = -(sum(Y(index1,:))-Y(index1,index1));
end


% SelfConductance definition
Y(8,8) = Y(8,8) + 1/(Rrby/2 + Rcov_i);
Y(9,9)= Y(9,9)+1/(Rh_slot/2 + Rcovo_slot);
Y(10,10)= Y(10,10)+1/(Rh_tooth/2 + Rcovo_tooth);

% Iteration calculation

% dT = zeros(1,10);
% dQ = zeros(1,10);
% if( FlagCon == 0)
% j = 1;
% T = zeros(1,10);
% end
% Z = zeros(10,1);
% Z(:,1) = 1;
% while(j<1600000)
%     
%     if j>1
%     Q = Z*T(j-1,:).*Y;
%     dQ(j,:) = P-sum(Q');
%     dT(j,:) = dQ(j,:)*delta./C;
%     else if j == 1
%             Q  = Z*[0 0 0 0 0 0 0 0 0 0].*Y;
%             dQ = P-sum(Q');
%             dT = dQ*delta./C;
%         end
%     end
%     if j > 1
%         T(j,:) = T(j-1,:)+dT(j,:);
%        if(abs(max(max(T)))>1)
%            break;
%        end
%         else if j == 1
%             T(j,:) = dT+T;
%         end
%     end
%     j=j+1;
%    if mod(j,5)==0
%        j
% 
%    end
% end





% dT = zeros(1,10);
% dQ = zeros(1,10);
% if( FlagCon == 0)
% j = 1;
% T = zeros(1,10);
% end
% Z = zeros(10,1);
% Z(:,1) = 1;
% while(j<1600000)
%     
%     if j>1
%     Q1 = Z*T(j-1,:).*Y;
%     dQ1 = P-sum(Q1');
%     dT1 = dQ1*delta./C/2;
%  
%     Q2 = Z*(T(j-1,:)+dT1).*Y;
%     dQ2 = P-sum(Q2');
%     dT2 = dQ2*delta./C/2;
%     
%     
%     else if j == 1
%             Q1  = Z*[0 0 0 0 0 0 0 0 0 0].*Y;
%             dQ1 = P-sum(Q1');
%             dT1 = dQ1*delta./C/2;
%             
%             Q2 = Z*([0 0 0 0 0 0 0 0 0 0]+dT1).*Y;
%             dQ2 = P-sum(Q2');
%             dT2 = dQ2*delta./C/2;
%         end
%     end
%     if j > 1
%         T(j,:) = T(j-1,:)+dT1+dT2;
%         else if j == 1
%             T(j,:) = dT1+dT2+T;
%         end
%     end
%     j=j+1;
%    if mod(j,1000)==0
%        j
% 
%    end
% end




dT = zeros(1,10);
dQ = zeros(1,10);
if( FlagCon == 0)
j = 1;
T = zeros(1000000,m);
end
Z = zeros(10,1);
Z(:,1) = 1;
while(j<1000001)
    
    if j>1
    Q1 = Z*T(j-1,:).*Y;
    dQ1 = P-sum(Q1');
    dT1 = dQ1*delta./C;
 
    Q2 = Z*(T(j-1,:)+dT1).*Y;
    dQ2 = P-sum(Q2');
    dT2 = dQ2*delta./C/2;
    
    Q3 = Z*(T(j-1,:)+dT1+dT2).*Y;
    dQ3 = P-sum(Q3');
    dT3 = dQ3*delta./C/2;

    Q4 = Z*(T(j-1,:)+dT1+dT2+dT3).*Y;
    dQ4 = P-sum(Q4');
    dT4 = dQ3*delta./C;
    
    else if j == 1
            Q1  = Z*[0 0 0 0 0 0 0 0 0 0].*Y;
            dQ1 = P-sum(Q1');
            dT1 = dQ1*delta./C;
            
            Q2 = Z*([0 0 0 0 0 0 0 0 0 0]+dT1).*Y;
            dQ2 = P-sum(Q2');
            dT2 = dQ2*delta./C/2;
            
            Q3 = Z*([0 0 0 0 0 0 0 0 0 0]+dT1+dT2).*Y;
            dQ3 = P-sum(Q3');
            dT3 = dQ3*delta./C/2;
            
            Q4 = Z*([0 0 0 0 0 0 0 0 0 0]+dT1+dT2+dT3).*Y;
            dQ4 = P-sum(Q4');
            dT4 = dQ3*delta./C;
        end
    end
    if j > 1
        T(j,:) = T(j-1,:)+dT1/6+dT2/3+dT3/3+dT4/6;
        else if j == 1
            T(j,:) = dT1/6+dT2/3+T(1,:)+dT3/3+dT4/6;
        end
    end
    j=j+1;
   if mod(j,5000)==0
       j

   end
end



C = [0.5 0 0.5];
time = (1:1:1000000)*delta/3600;
plot(time,T(:,1),'-k','LineWidth',1.5);
hold on
plot(time,T(:,2),'--r','LineWidth',1.5);
plot(time,T(:,3),':','LineWidth',1.5,'Color',C);
plot(time,T(:,4),'-.r','LineWidth',1.5);
plot(time,T(:,5),':k','LineWidth',1.5);
plot(time,T(:,6),'--','LineWidth',1.5,'Color',C);
plot(time,T(:,7),':g','LineWidth',1.5);
plot(time,T(:,8),'-','LineWidth',1.5,'Color',C);
plot(time,T(:,9),'-r','LineWidth',1.5);
plot(time,T(:,10),'-b','LineWidth',1.5);
legend('Winding','Tooth','Back iron slot','Back iron tooth','Airgap','Rotor slot','Rotor tooth','Rotor back iron','Housing slot','Housing tooth');
ylabel('Temperature rise(K)');
xlabel('Time(Hour)');
set(gca, 'Fontname', 'Times New Roman','FontSize',24);
%legend('6^{th} harmonic open-circuit','12^{th} harmonic open-circuit','Half pkk indcued voltage')
%set(gca,'XTick',[0:Q/6:Q]) %改变x轴坐标间隔显示 这里间隔为2
%set(gca,'YTick',[0:0.2:1]) %改变x轴坐标间隔显示 这里间隔为2
%set(gca,'YTick',[0:10:60]) %改变x轴坐标间隔显示 这里间隔为2
set(gcf,'Position',[0/0.277 0/0.277 200/0.277 160/0.277]);
set(gca,'linewidth',1.5);

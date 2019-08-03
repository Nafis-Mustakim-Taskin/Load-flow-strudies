clear all
close all
clc
n=6; % total number of buses
m=2; %number of pv buses
a=m+1;

Z1 = zeros (n,n);
Z2 = zeros (n,n);
Z = zeros (n,n);
Y=  zeros (n,n);
Ybus= zeros (n,n);
%Impedence value input 
for x=1:n
    for y=1:n
        if y>x 
            prompt = ' Give values '; 
            sprintf('Impedence from bus %d to bus %d',x,y)
            Z1(x,y) = input(prompt);
        end
    end
end
for x=1:n
    for y=1:n
       if x>y
        Z2(x,y) = Z1(y,x);
       end
       
       end
end
Z = Z1 + Z2 ; % Final impedence matrix
% Forming Y Bus 
for x=1:n
    for y=1:n 
        if x~=y && Z(x,y)~= 0
            Y(x,y) = Z(x,y).^-1;
        end
    end
end
for x=1:1:n
    for y=1:1:n
        if x==y 
            for k=1:1:n
                if k~=x
                    Y(x,y)=Y(x,y)-Y(x,k);
                end
            end
        end
    end
end
Ybus = Y .*i;

P1=0;
Q1=0;
P2=0;
P3=0;
Q2=0;
Q3=0;
P4=0;
P5=0;
Q4=0;
Q5=0;
P6=0;
Q6=0;

Pspec=[P1 .8 .5 .1 .08 .1]; %Active power should be specified here
Qspec=[ Q1 Q2 Q3 .04 .03 .03]; % Reactive power for PQ Buses should be specified for reactive buses Q index should be shown
Voltagearrayold=[1 1 1 1 1 1]; % The first entry is for slack bus. Then voltages for PV buses and PQ buses ( initial assumptions) should be specified respectively.
Deltaarrayold=[0 0 0 0 0 0]; % Estimated Angles for the buses should be specified here.
angleslackbus=Deltaarrayold(1);
Ybusmagnitude=abs(Ybus);
Ybustheta=angle(Ybus);

ni=10; % Number of iterations

for j=1:ni
%Initialisation    
P1=0;
Q1=0;
P2=0;
P3=0;
Q2=0;
Q3=0;
P4=0;
P5=0;
P6=0;
Q4=0;
Q5=0;
Q6=0;

% P and Q calculation   
for i=1:n
         
         P1=P1+(Voltagearrayold(1)*Voltagearrayold(i)*Ybusmagnitude(1,i))*(cos(-Ybustheta(1,i)-Deltaarrayold(i)+Deltaarrayold(1)));
        
         Q1=(Q1+(Voltagearrayold(1)*Voltagearrayold(i)*Ybusmagnitude(1,i))*(sin(-Ybustheta(1,i)-Deltaarrayold(i)+Deltaarrayold(1))));
         P2=P2+(Voltagearrayold(2)*Voltagearrayold(i)*Ybusmagnitude(2,i))*(cos(-Ybustheta(2,i)-Deltaarrayold(i)+Deltaarrayold(2)));
         Q2=(Q2+(Voltagearrayold(2)*Voltagearrayold(i)*Ybusmagnitude(2,i))*(sin(-Ybustheta(2,i)-Deltaarrayold(i)+Deltaarrayold(2))));
         P3=P3+(Voltagearrayold(3)*Voltagearrayold(i)*Ybusmagnitude(3,i))*(cos(-Ybustheta(3,i)-Deltaarrayold(i)+Deltaarrayold(3)));
         Q3=Q3+(Voltagearrayold(3)*Voltagearrayold(i)*Ybusmagnitude(3,i))*(sin(-Ybustheta(3,i)-Deltaarrayold(i)+Deltaarrayold(3)));
         P4=P4+(Voltagearrayold(4)*Voltagearrayold(i)*Ybusmagnitude(4,i))*(cos(-Ybustheta(4,i)-Deltaarrayold(i)+Deltaarrayold(4)));
         P5=P5+(Voltagearrayold(5)*Voltagearrayold(i)*Ybusmagnitude(5,i))*(cos(-Ybustheta(5,i)-Deltaarrayold(i)+Deltaarrayold(5)));
         P6=P6+(Voltagearrayold(6)*Voltagearrayold(i)*Ybusmagnitude(6,i))*(cos(-Ybustheta(6,i)-Deltaarrayold(i)+Deltaarrayold(6)));
         Q4=Q4+(Voltagearrayold(4)*Voltagearrayold(i)*Ybusmagnitude(4,i))*(sin(-Ybustheta(4,i)-Deltaarrayold(i)+Deltaarrayold(4)));
         Q5=Q5+(Voltagearrayold(5)*Voltagearrayold(i)*Ybusmagnitude(5,i))*(sin(-Ybustheta(5,i)-Deltaarrayold(i)+Deltaarrayold(5)));
         Q6=Q6+(Voltagearrayold(6)*Voltagearrayold(i)*Ybusmagnitude(6,i))*(sin(-Ybustheta(6,i)-Deltaarrayold(i)+Deltaarrayold(6)));
end
%disp(Q2);
%disp(Q3);
%disp(Q4);
%disp(Q5);
     % Calculation of the mismatches
     MISMATCHP2=(Pspec(2))-P2;
     MISMATCHP3=(Pspec(3))-P3;
     MISMATCHP4=(Pspec(4))-P4;
     MISMATCHP5=(Pspec(5))-P5;
     MISMATCHP6=(Pspec(6))-P6;
     MISMATCHQ2=(Qspec(2))-Q2;
 %    disp('Test');
  %   disp(MISMATCHQ2);
     MISMATCHQ3=(Qspec(3))-Q3;
    % disp(MISMATCHQ3);
     MISMATCHQ4=(Qspec(4))-Q4;
     %disp(MISMATCHQ4);
     MISMATCHQ5=(Qspec(5))-Q5;
   %  disp(MISMATCHQ5);
     MISMATCHQ6=(Qspec(6))-Q6;
     MISMATCHMATRIX=[MISMATCHP2 MISMATCHP3 MISMATCHP4 MISMATCHP5 MISMATCHP6  MISMATCHQ2 MISMATCHQ3 MISMATCHQ4 MISMATCHQ5 MISMATCHQ6 ];
     %disp(n);
     %disp(a);
     MISMATCHMATRIX(n:n+m-1)=[];
     
     disp(MISMATCHMATRIX);
     % Formation of the mismatch matrix
     MISMATCH=MISMATCHMATRIX.';
     % Formation of the Jacobian
     %Initialization
 J1=zeros((n-1),(n-1));
J2=zeros((n-1),(n-a));
J3=zeros((n-a),(n-1));
J4=zeros((n-a),(n-a));
     for x=2:n
         for y=2:n
             if x==y
             for k=1:n
                 if k==x continue
                 end
                 
             J1(x-1,y-1)=J1(x-1,y-1)+(-Voltagearrayold(x)*Voltagearrayold(k)*Ybusmagnitude(x,k)*sin(Deltaarrayold(x)-Deltaarrayold(k)-Ybustheta(x,k)));
             end
             
             else J1(x-1,y-1)=Voltagearrayold(x)*Voltagearrayold(y)*Ybusmagnitude(x,y)*sin(Deltaarrayold(x)-Deltaarrayold(y)-Ybustheta(x,y));
         end
     end
     
     end
     for x=2:n
         for y=a+1:n
             if x==y
             for k=1:n
                 if k==x continue
                 end
                 J2(x-1,y-a)=J2(x-1,y-a)+(Voltagearrayold(k)*Ybusmagnitude(x,k)*cos(Deltaarrayold(x)-Deltaarrayold(k)-Ybustheta(x,k)));
             end 
             J2(x-1,y-a)=2*Voltagearrayold(x)*real(Ybus(x,x))+J2(x-1,y-a);
             else J2(x-1,y-a)=Voltagearrayold(x)*Ybusmagnitude(x,y)*cos(Deltaarrayold(x)-Deltaarrayold(y)-Ybustheta(x,y));
             end
         end
     end
     for x=a+1:n
         for y=2:n
             if x==y
             for k=1:n
                 if k==x continue
                 end
                 J3(x-a,y-1)=J3(x-a,y-1)+(Voltagearrayold(x)*Voltagearrayold(k)*Ybusmagnitude(x,k)*cos(Deltaarrayold(x)-Deltaarrayold(k)-Ybustheta(x,k)));
             end
             else J3(x-a,y-1)=-Voltagearrayold(x)*Voltagearrayold(y)*Ybusmagnitude(x,y)*cos(Deltaarrayold(x)-Deltaarrayold(y)-Ybustheta(x,y));
             end 
         end
         
     end
     for x=a+1:n
         for y=a+1:n
             if x==y
             for k=1:n
                 if k==x continue
                 end
                 
                 J4(x-a,y-a)=J4(x-a,y-a)+(Voltagearrayold(k)*Ybusmagnitude(x,k)*sin(Deltaarrayold(x)-Deltaarrayold(k)-Ybustheta(x,k)));
             end
             J4(x-a,y-a)=-2*Voltagearrayold(x)*imag(Ybus(x,x))+ J4(x-a,y-a);
             else J4(x-a,y-a)=Voltagearrayold(x)*Ybusmagnitude(x,y)*sin(Deltaarrayold(x)-Deltaarrayold(y)-Ybustheta(x,y));
             end
         end
     end
     % J corresponds to the Jacobian matrix.
     J=[J1 J2;
         J3 J4];
     %disp(J);
     % Calculation of the correction vector
     CorrectionVector=J\MISMATCH;
     %disp(CorrectionVector);
     
     Deltaarrayintermediate=(Deltaarrayold(2:n)).'+CorrectionVector(1:n-1);
     Deltaarraynew=vertcat(angleslackbus.',Deltaarrayintermediate);
     % Output of angles
     disp('All the angles of the buses');
     disp(Deltaarraynew);
          Voltagearrayintermediate=CorrectionVector(n:2*n-m-2)+(Voltagearrayold(a+1:n)).';
          V1=Voltagearrayintermediate.';
         % disp(Voltagearrayintermediate);
     Voltagearraynew=vertcat((Voltagearrayold(1:a)).',(V1)');
     %Output of voltages
     disp ('All the voltages of the buses');
     disp(Voltagearraynew);
     Voltagearrayold=Voltagearraynew.';
     Deltaarrayold=Deltaarraynew.';
      for i=1:n
     P1=P1+(Voltagearrayold(1)*Voltagearrayold(i)*Ybusmagnitude(1,i))*(cos(-Ybustheta(1,i)-Deltaarrayold(i)+Deltaarrayold(1)));
      
         Q1=(Q1+(Voltagearrayold(1)*Voltagearrayold(i)*Ybusmagnitude(1,i))*(sin(-Ybustheta(1,i)-Deltaarrayold(i)+Deltaarrayold(1))));
         P2=P2+(Voltagearrayold(2)*Voltagearrayold(i)*Ybusmagnitude(2,i))*(cos(-Ybustheta(2,i)-Deltaarrayold(i)+Deltaarrayold(2)));
         Q2=(Q2+(Voltagearrayold(2)*Voltagearrayold(i)*Ybusmagnitude(2,i))*(sin(-Ybustheta(2,i)-Deltaarrayold(i)+Deltaarrayold(2))));
         P3=P3+(Voltagearrayold(3)*Voltagearrayold(i)*Ybusmagnitude(3,i))*(cos(-Ybustheta(3,i)-Deltaarrayold(i)+Deltaarrayold(3)));
         Q3=Q3+(Voltagearrayold(3)*Voltagearrayold(i)*Ybusmagnitude(3,i))*(sin(-Ybustheta(3,i)-Deltaarrayold(i)+Deltaarrayold(3)));
         P4=P4+(Voltagearrayold(4)*Voltagearrayold(i)*Ybusmagnitude(4,i))*(cos(-Ybustheta(4,i)-Deltaarrayold(i)+Deltaarrayold(4)));
         P5=P5+(Voltagearrayold(5)*Voltagearrayold(i)*Ybusmagnitude(5,i))*(cos(-Ybustheta(5,i)-Deltaarrayold(i)+Deltaarrayold(5)));
         P6=P6+(Voltagearrayold(6)*Voltagearrayold(i)*Ybusmagnitude(6,i))*(cos(-Ybustheta(6,i)-Deltaarrayold(i)+Deltaarrayold(6)));
         Q4=Q4+(Voltagearrayold(4)*Voltagearrayold(i)*Ybusmagnitude(4,i))*(sin(-Ybustheta(4,i)-Deltaarrayold(i)+Deltaarrayold(4)));
         Q5=Q5+(Voltagearrayold(5)*Voltagearrayold(i)*Ybusmagnitude(5,i))*(sin(-Ybustheta(5,i)-Deltaarrayold(i)+Deltaarrayold(5)));
         Q6=Q6+(Voltagearrayold(6)*Voltagearrayold(i)*Ybusmagnitude(6,i))*(sin(-Ybustheta(6,i)-Deltaarrayold(i)+Deltaarrayold(6)));
     end
ActivePowerMatrix=[P1 P2 P3 P4 P5 P6]; 
disp('Active Power Values');
disp(ActivePowerMatrix);
ReactivePowerMatrix=[Q1 Q2 Q3 Q4 Q5 Q6]; 
disp('Reactive Power Values');
disp(ReactivePowerMatrix);
                 end
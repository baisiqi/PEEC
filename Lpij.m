function Lp = Lpij( w1,w2,Dx,l1,l2,Dy,Dz )
%Lp is nH
% partial inductance calculation
SMALL_NO=1e-10;
a(1)=Dx-l1/2-l2/2;
a(2)=Dx+l1/2-l2/2;
a(3)=Dx+l1/2+l2/2;
a(4)=Dx-l1/2+l2/2;
b(1)=Dy-w1/2-w2/2;
b(2)=Dy+w1/2-w2/2;
b(3)=Dy+w1/2+w2/2;
b(4)=Dy-w1/2+w2/2;

Lp=0;

for i=1:4
    for j=1:4
        p=(a(i)^2+b(j)^2+Dz^2)^0.5+SMALL_NO;   
        A=(b(j)^2-Dz^2)/2*a(i)*log(a(i)+p);
        B=(a(i)^2-Dz^2)/2*b(j)*log(b(j)+p);
        C=-(b(j)^2-2*Dz^2+a(i)^2)*p/6;
        D=-b(j)*Dz*a(i)*atan(a(i)*b(j)/p/(Dz+SMALL_NO));
        Lp=Lp+(A+B+C+D)*(-1)^(i+j);
    end
end

Lp=Lp*0.1/w1/w2;
end

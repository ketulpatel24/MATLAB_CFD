clc
clear
clear all:
X=1;Y=1;Z=1;
Mx=32;Ny=32;Lz=32;
A=[log2(Mx),log2(Ny),log2(Lz)];
level=min(A);
%level=3;
m=zeros(level);n=zeros(level);l=zeros(level);
dx=zeros(level);dy=zeros(level);dz=zeros(level);
presweep=12;postsweep=12;
%data=fopen('data.txt','w');

%calculate m,n,l values at each level
for i=1:level
    m(i)=Mx/power(2,i-1);n(i)=Ny/power(2,i-1);l(i)=Lz/power(2,i-1);
end
%calculate dx,dy,dz values at each level
for i=1:level
    dx(i)=X/m(i);dy(i)=Y/n(i);dz(i)=Z/l(i);
end
%define phi for each level
for N=1:level
    for i=1:m(N)+2
        for j=1:n(N)+2
            for k=1:l(N)+2
                phi(N,i,j,k)=0;
                phi0(N,i,j,k)=0;
                %phip(N,i,j,k)=0;
                %res(N,i,j,k)=0;
            end
        end
    end
    
end

%define coefficients for each level
for N=1:level
    for i=1:m(N)+2
        for j=1:n(N)+2
            for k=1:l(N)+2
                ae(N,i,j,k)=1/power(dx(N),2);aw(N,i,j,k)=1/power(dx(N),2);
                an(N,i,j,k)=1/power(dy(N),2);as(N,i,j,k)=1/power(dy(N),2);
                at(N,i,j,k)=1/power(dz(N),2);ab(N,i,j,k)=1/power(dz(N),2);
                ap(N,i,j,k)=ae(N,i,j,k)+aw(N,i,j,k)+an(N,i,j,k)+as(N,i,j,k)+at(N,i,j,k)+ab(N,i,j,k);
                
            end
        end
    end
    
    
end

%calculate x,y,z values for each level
for N=1:level
    x(N,1)=-dx(N)/2;x(N,m(N)+2)=1+(dx(N)/2);  %define Ghost Points
    y(N,1)=-dy(N)/2;y(N,n(N)+2)=1+(dy(N)/2);
    z(N,1)=-dz(N)/2;z(N,l(N)+2)=1+(dz(N)/2);
    for i=2:m(N)+1
        x(N,i)=x(N,i-1)+dx(N);
    end
    for j=2:n(N)+1
        y(N,j)=y(N,j-1)+dy(N);
    end
    for k=2:l(N)+1
        z(N,k)=z(N,k-1)+dz(N);
    end
end

%calculate RHS term for first level
for k=1:l(1)+2
    for j=1:n(1)+2
        for i=1:m(1)+2
                b(i,j,k)=50000*exp(-50*(power(1-x(1,i),2)+power(z(1,k),2)))*(100*(power(1-x(1,i),2)+power(z(1,k),2))-2);
        end
    end
end

L2norm=1.0;count=0;
while(count<100)
    %L2norm>0.00001
    %storing phi values in phi0 at level 1
    for k=1:l(1)+2
        for j=1:n(1)+2
            for i=1:m(1)+2
                phi0(1,i,j,k)=phi(1,i,j,k);
            end
        end
    end
    
    %gauss seidel loop for first level
    for t=1:presweep
    %calculation of phi
        for k=2:l(1)+1
            for j=2:n(1)+1
                for i=2:m(1)+1                 
                    phi(1,i,j,k)=(ae(1,i,j,k)*phi(1,i+1,j,k)+aw(1,i,j,k)*phi(1,i-1,j,k)+an(1,i,j,k)*phi(1,i,j+1,k) ...
                        +as(1,i,j,k)*phi(1,i,j-1,k)+at(1,i,j,k)*phi(1,i,j,k+1)+ab(1,i,j,k)*phi(1,i,j,k-1)-b(i,j,k))/ap(1,i,j,k);

                end
            end
        end
    %Update boundary and Periodic condition
        for k=1:l(1)+2
            for j=1:n(1)+2
                for i=1:m(1)+2
                    phi(1,m(1)+2,j,k)=2*(100*(1-z(1,k))+500*exp(-50*power(z(1,k),2)))-phi(1,m(1)+1,j,k);
                    phi(1,1,j,k)=2*(500*exp(-50*(1+power(z(1,k),2))))-phi(1,2,j,k);
                    phi(1,i,j,1)=2*(100*x(1,i)+500*exp(-50*power(1-x(1,i),2)))-phi(1,i,j,2);
                    phi(1,i,j,l(1)+2)=2*(500*exp(-50*(power(1-x(1,i),2)+1)))-phi(1,i,j,l(1)+1);
                
                    phi(1,i,n(1)+2,k)=phi(1,i,2,k);
                    phi(1,i,1,k)=phi(1,i,n(1)+1,k);
                end
            end
        end
    end
    
    for k=1:l(1)+2
        for i=1:m(1)+2
            v1(k,i)=phi(1,i,(n(1)+2)/2,k);
        
        %fprintf(data,'%6.2f %6.2f %12.8f\r\n',x(1,i),z(1,k),u);
        %v(k,i)=500*exp(-50*(((1-x(1,i))^2)+z(1,k)*z(1,k)))+100*x(1,i)*(1-z(1,k));
        end
    end
    %calculation of residual at 1st level
    for k=2:l(1)+1
        for j=2:n(1)+1
            for i=2:m(1)+1
                res(1,i,j,k)=b(i,j,k)-(ae(1,i,j,k)*phi(1,i+1,j,k)+aw(1,i,j,k)*phi(1,i-1,j,k)+an(1,i,j,k)*phi(1,i,j+1,k)+ ...
                    as(1,i,j,k)*phi(1,i,j-1,k)+at(1,i,j,k)*phi(1,i,j,k+1)+ab(1,i,j,k)*phi(1,i,j,k-1)-ap(1,i,j,k)*phi(1,i,j,k));
            end
        end
    end
    %update residual at boundary
    for k=1:l(1)+2
        for j=1:n(1)+2
            for i=1:m(1)+2
                res(1,m(1)+2,j,k)=-res(1,m(1)+1,j,k);
                res(1,1,j,k)=-res(1,2,j,k);
                res(1,i,j,1)=-res(1,i,j,2);
                res(1,i,j,l(1)+2)=-res(1,i,j,l(1)+1);
                res(1,i,1,k)=-res(1,i,2,k);
                res(1,i,n(1)+2,k)=-res(1,i,n(1)+1,k);
            end
        end
    end
    
    for N=2:level
        
        %restrict residual of current level to next level RHS
        for k=2:l(N)+1
            for j=2:n(N)+1
                for i=2:m(N)+1
                    p=2*i-2;q=2*j-2;r=2*k-2;
                    rhs(N,i,j,k)=0.125*(res(N-1,p,q,r)+res(N-1,p+1,q,r)+res(N-1,p,q+1,r)+res(N-1,p+1,q+1,r)+res(N-1,p,q,r+1)+ ...
                        res(N-1,p+1,q,r+1)+res(N-1,p,q+1,r+1)+res(N-1,p+1,q+1,r+1));
                end
            end
        end
        
        %update RHS values at boundary
        for k=1:l(N)+2
            for j=1:n(N)+2
                for i=1:m(N)+2
                    rhs(N,m(N)+2,j,k)=-rhs(N,m(N)+1,j,k);
                    rhs(N,1,j,k)=-rhs(N,2,j,k);
                    rhs(N,i,j,1)=-rhs(N,i,j,2);
                    rhs(N,i,j,l(N)+2)=-rhs(N,i,j,l(N)+1);
                    rhs(N,i,1,k)=-rhs(N,i,2,k);
                    rhs(N,i,n(N)+2,k)=-rhs(N,i,n(N)+1,k);
                end
            end
        end
        
        %Initialize phi as error for next level equal to zero
        for k=1:l(N)+2
            for j=1:n(N)+2
                for i=1:m(N)+2
                    phi(N,i,j,k)=0;
                end
            end
        end
        
        
        %GS loop for each level except 1st level
        for t=1:presweep
        %calculation of phi at each level except 1st level
            for k=2:l(N)+1
                for j=2:n(N)+1
                    for i=2:m(N)+1
                        phi(N,i,j,k)=(ae(N,i,j,k)*phi(N,i+1,j,k)+aw(N,i,j,k)*phi(N,i-1,j,k)+an(N,i,j,k)*phi(N,i,j+1,k)+ ...
                            as(N,i,j,k)*phi(N,i,j-1,k)+at(N,i,j,k)*phi(N,i,j,k+1)+ab(N,i,j,k)*phi(N,i,j,k-1)-rhs(N,i,j,k))/ap(N,i,j,k);
                    end
                end
            end
            
            %update boundary condition at each level
            for k=1:l(N)+2
                for j=1:n(N)+2
                    for i=1:m(N)+2
                        phi(N,m(N)+2,j,k)=-phi(N,m(N)+1,j,k);
                        phi(N,1,j,k)=-phi(N,2,j,k);
                        phi(N,i,j,1)=-phi(N,i,j,2);
                        phi(N,i,j,l(N)+2)=-phi(N,i,j,l(N)+1);
                
                        phi(N,i,n(N)+2,k)=-phi(N,i,n(N)+1,k);
                        phi(N,i,1,k)=-phi(N,i,2,k);
                    end
                end
            end
        end
        %calculation of residual at each level
        for k=2:l(N)+1
            for j=2:n(N)+1
                for i=2:m(N)+2
                    res(N,i,j,k)=rhs(N,i,j,k)-(ae(N,i,j,k)*phi(N,i+1,j,k)+aw(N,i,j,k)*phi(N,i-1,j,k)+an(N,i,j,k)*phi(N,i,j+1,k)+ ...
                    as(N,i,j,k)*phi(N,i,j-1,k)+at(N,i,j,k)*phi(N,i,j,k+1)+ab(N,i,j,k)*phi(N,i,j,k-1)-ap(N,i,j,k)*phi(N,i,j,k));
                end
            end
        end
        
        %update residual at boundary
        for k=1:l(N)+2
            for j=1:n(N)+2
                for i=1:m(N)+2
                    res(N,m(N)+2,j,k)=-res(N,m(N)+1,j,k);
                    res(N,1,j,k)=-res(N,2,j,k);
                    res(N,i,j,1)=-res(N,i,j,2);
                    res(N,i,j,l(N)+2)=-res(N,i,j,l(N)+1);
                    res(N,i,1,k)=-res(N,i,2,k);
                    res(N,i,n(N)+2,k)=-res(N,i,n(N)+1,k);
                end
            end
        end
    end
    
    
    for k=2:l(2)+1
        for i=1:m(2)+1
                v2(k,i)=phi(2,i,(n(2)+2)/2,k);
        end
    end
    for k=2:l(3)+1
        for i=2:m(3)+1
            v3(k,i)=phi(3,i,(n(3)+2)/2,k);
        end
    end
    for k=2:l(4)+1
        for i=2:m(4)+1
            v4(k,i)=phi(4,i,(n(4)+2)/2,k);
        end
    end
    for k=2:l(5)+1
        for i=1:m(5)+2
            v5(k,i)=phi(5,i,(n(4)+2)/2,k);
        end
    end
    
    for N=level:2
        %prolongate phi into phi(N)(N-1)
         for k=2:l(N)+1
             for j=2:n(N)+1
                 for i=2:m(N)+1
                     p=2*i-1;q=2*j-1;r=2*k-1;
                        phip(N,p,q,r)=(1/64)*(27*phi(N,i,j,k)+9*phi(N,i+1,j,k)+9*phi(N,i,j+1,k)+9*phi(N,i,j,k+1)+3*phi(N,i+1,j+1,k)+ ...
                            3*phi(N,i+1,j,k+1)+3*phi(N,i,j+1,k+1)+phi(N,i+1,j+1,k+1));
                        phip(N,p+1,q,r)=(1/64)*(9*phi(N,i,j,k)+27*phi(N,i+1,j,k)+3*phi(N,i,j+1,k)+3*phi(N,i,j,k+1)+9*phi(N,i+1,j+1,k)+ ...
                            9*phi(N,i+1,j,k+1)+phi(N,i,j+1,k+1)+3*phi(N,i+1,j+1,k+1));
                        phip(N,p,q+1,r)=(1/64)*(9*phi(N,i,j,k)+3*phi(N,i+1,j,k)+27*phi(N,i,j+1,k)+3*phi(N,i,j,k+1)+9*phi(N,i+1,j+1,k)+ ...
                            phi(N,i+1,j,k+1)+9*phi(N,i,j+1,k+1)+3*phi(N,i+1,j+1,k+1));
                        phip(N,p,q,r+1)=(1/64)*(9*phi(N,i,j,k)+3*phi(N,i+1,j,k)+3*phi(N,i,j+1,k)+27*phi(N,i,j,k+1)+phi(N,i+1,j+1,k)+ ...
                            9*phi(N,i+1,j,k+1)+9*phi(N,i,j+1,k+1)+3*phi(N,i+1,j+1,k+1));
                        phip(N,p+1,q,r+1)=(1/64)*(3*phi(N,i,j,k)+9*phi(N,i+1,j,k)+phi(N,i,j+1,k)+9*phi(N,i,j,k+1)+3*phi(N,i+1,j+1,k)+ ...
                            27*phi(N,i+1,j,k+1)+3*phi(N,i,j+1,k+1)+9*phi(N,i+1,j+1,k+1));
                        phip(N,p,q+1,r+1)=(1/64)*(3*phi(N,i,j,k)+phi(N,i+1,j,k)+9*phi(N,i,j+1,k)+9*phi(N,i,j,k+1)+9*phi(N,i+1,j+1,k)+ ...
                            3*phi(N,i+1,j,k+1)+27*phi(N,i,j+1,k+1)+3*phi(N,i+1,j+1,k+1));
                        phip(N,p+1,q+1,r)=(1/64)*(3*phi(N,i,j,k)+9*phi(N,i+1,j,k)+9*phi(N,i,j+1,k)+phi(N,i,j,k+1)+27*phi(N,i+1,j+1,k)+ ...
                            3*phi(N,i+1,j,k+1)+3*phi(N,i,j+1,k+1)+9*phi(N,i+1,j+1,k+1));
                        phip(N,p+1,q+1,r+1)=(1/64)*(phi(N,i,j,k)+3*phi(N,i+1,j,k)+3*phi(N,i,j+1,k)+3*phi(N,i,j,k+1)+9*phi(N,i+1,j+1,k)+ ...
                            9*phi(N,i+1,j,k+1)+9*phi(N,i,j+1,k+1)+27*phi(N,i+1,j+1,k+1));
                 end
             end
         end
         
         for k=1:l(N)+2
            for j=1:n(N)+2
                for i=1:m(N)+2
                    phip(N,m(N)+2,j,k)=-phip(N,m(N)+1,j,k);
                    phip(N,1,j,k)=-phip(N,2,j,k);
                    phip(N,i,j,1)=-phip(N,i,j,2);
                    phip(N,i,j,l(N)+2)=-phip(N,i,j,l(N)+1);
                    phip(N,i,1,k)=-phip(N,i,2,k);
                    phip(N,i,n(N)+2,k)=-phip(N,i,n(N)+1,k);
                end
            end
        end
         %update value of phi1
        for k=1:l(N-1)+2
            for j=1:n(N-1)+2
                for i=1:m(N-1)+2
                    phi(N-1,i,j,k)=phi(N-1,i,j,k)+phip(N,i,j,k);
                end
            end
        end
        
        %GS loop for previous level
        for t=1:postsweep
        %calculation of phi at previous level
            for k=2:l(N-1)+1
                for j=2:n(N-1)+1
                    for i=2:m(N-1)+1
                        phi(N-1,i,j,k)=(ae(N-1,i,j,k)*phi(N-1,i+1,j,k)+aw(N-1,i,j,k)*phi(N-1,i-1,j,k)+an(N-1,i,j,k)*phi(N-1,i,j+1,k)+ ...
                            as(N-1,i,j,k)*phi(N-1,i,j-1,k)+at(N-1,i,j,k)*phi(N-1,i,j,k+1)+ab(N-1,i,j,k)*phi(N-1,i,j,k-1)-rhs(N-1,i,j,k))/ap(N-1,i,j,k);
                    end
                end
            end
    
            %update boundary condition
            for k=1:l(N-1)+2
                for j=1:n(N-1)+2
                    for i=1:m(N-1)+2
                        phi(N-1,m(N-1)+2,j,k)=-phi(N-1,m(N-1)+1,j,k);
                        phi(N-1,1,j,k)=-phi(N-1,2,j,k);
                        phi(N-1,i,j,1)=-phi(N-1,i,j,2);
                        phi(N-1,i,j,l(N-1)+2)=-phi(N-1,i,j,l(N-1)+1);
                
                        phi(N-1,i,n(N-1)+2,k)=-phi(N-1,i,n(N-1)+1,k);
                        phi(N-1,i,1,k)=-phi(N-1,i,2,k);
                    end
                end
            end
        end

    end
    
    %gauss seidel loop for first level
    for t=1:postsweep
        %calculation of phi
        for k=2:l(1)+1
            for j=2:n(1)+1
                for i=2:m(1)+1                 
                    phi(1,i,j,k)=(ae(1,i,j,k)*phi(1,i+1,j,k)+aw(1,i,j,k)*phi(1,i-1,j,k)+an(1,i,j,k)*phi(1,i,j+1,k) ...
                        +as(1,i,j,k)*phi(1,i,j-1,k)+at(1,i,j,k)*phi(1,i,j,k+1)+ab(1,i,j,k)*phi(1,i,j,k-1)-b(i,j,k))/ap(1,i,j,k);

                end
            end
        end
        %Update boundary and Periodic condition
        for k=1:l(1)+2
            for j=1:n(1)+2
                for i=1:m(1)+2
                    phi(1,m(1)+2,j,k)=2*(100*(1-z(1,k))+500*exp(-50*power(z(1,k),2)))-phi(1,m(1)+1,j,k);
                    phi(1,1,j,k)=2*(500*exp(-50*(1+power(z(1,k),2))))-phi(1,2,j,k);
                    phi(1,i,j,1)=2*(100*x(1,i)+500*exp(-50*power(1-x(1,i),2)))-phi(1,i,j,2);
                    phi(1,i,j,l(1)+2)=2*(500*exp(-50*(power(1-x(1,i),2)+1)))-phi(1,i,j,l(1)+1);
                
                    phi(1,i,n(1)+2,k)=phi(1,i,2,k);
                    phi(1,i,1,k)=phi(1,i,n(1)+1,k);
                end
            end
        end
    end
    for k=1:l(1)+2
        for i=1:m(1)+2
            u1(k,i)=phi(1,i,(n(1)+2)/2,k);
        
        %fprintf(data,'%6.2f %6.2f %12.8f\r\n',x(1,i),z(1,k),u);
        %v(k,i)=500*exp(-50*(((1-x(1,i))^2)+z(1,k)*z(1,k)))+100*x(1,i)*(1-z(1,k));
        end
    end
    for k=2:l(2)+1
        for i=2:m(2)+1
            u2(k,i)=phi(2,i,(n(2)+2)/2,k);
        end
    end
for k=2:l(3)+1
    for i=2:m(3)+1
        u3(k,i)=phi(3,i,(n(3)+2)/2,k);
    end
end
for k=2:l(4)+1
    for i=2:m(4)+1
        u4(k,i)=phi(4,i,(n(4)+2)/2,k);
    end
end
    
    sum=0;
    %calculate sum of square of error
    for k=1:l(1)+2
        for j=1:n(1)+2
            for i=1:m(1)+2
                sum=sum+power(phi(1,i,j,k)-phi0(1,i,j,k),2);
            end
        end
    end
    L2norm=sqrt(sum);
   
    count=count+1;

    disp(count)
    disp(L2norm)
end

    for k=1:l(1)+2
        for i=1:m(1)+2
            u(k,i)=phi(1,i,(n(1)+2)/2,k);
        
        %fprintf(data,'%6.2f %6.2f %12.8f\r\n',x(1,i),z(1,k),u);
        v(k,i)=500*exp(-50*(((1-x(1,i))^2)+z(1,k)*z(1,k)))+100*x(1,i)*(1-z(1,k));
        end
    end

%fclose(data);
g=[50,100,150,200,250,300,350,400,450,500];
% figure(1)
% contourf(v1,10)
% c.LineWidth = 2;
% title('level-1 after 100 V cycle ');
% figure(2)
% contourf(v2,10)
% c.LineWidth = 2;
% title('level-2 (restriction) after 100 V cycle ');
% figure(3)
% contourf(v3,10)
% c.LineWidth = 2;
% title('level-3 (restriction) after 100 V cycle ');
% figure(4)
% contourf(v4,10)
% c.LineWidth = 2;
% title('level-4 (restriction) after 100 V cycle ');
% figure(5)
% contourf(v5,10)
% c.LineWidth = 2;
% title('level-5 (restriction) after 100 V cycle ')
% figure(6)
% contourf(u4,10)
% c.LineWidth = 2;
% title('level-4 (prolongation) after 100 V cycle ');
% figure(7)
% contourf(u3,10)
% c.LineWidth = 2;
% title('level-3 (prolongation) after 100 V cycle ');
% figure(8)
% contourf(u2,10)
% c.LineWidth = 2;
% title('level-2 (prolongation) after 100 V cycle ');
figure(9)
contourf(u,10)
c.LineWidth = 2;
title('Contour plot for 32*32*32')
figure(10)
contourf(v,10)
c.LineWidth = 2;
title('Contour plot for Analytical solution')
% figure(11)
% contourf(u,10)
% c.LineWidth = 2;
% title('level-1 (prolongation) after 100 V cycle ');

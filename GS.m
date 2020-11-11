clc
clear

%define x,y,z length
Lx = 1;
Ly = 1;
Lz = 1;

%define number of grid in each direction
Nx = 32;  Ny = 32;  Nz = 32;

%calculate dx, dy,dz
dx = Lx/Nx; dy = Ly/Ny; dz = Lz/Nz;

%calculate x,y and z
x(1) = -dx/2; y(1) = -dy/2; z(1) = -dz/2;
x(Nx+2) = Lx + dx/2; y(Ny+2) = Ly + dy/2; z(Nz+2) = Lz + dz/2;

for i=2:Nx+1
    x(i) = x(i-1) + dx;
end

for j=2:Ny+1
    y(j) = y(j-1) + dy;
end

for k=2:Nz+1
    z(k) = z(k-1) + dz;
end

%calculate coefficients
for k=2:Nz+1
    for j=2:Ny+1
        for i=2:Nx+1
            ae(i,j,k) = 1/power(dx,2); aw(i,j,k) = 1/power(dx,2);
            an(i,j,k) = 1/power(dy,2); as(i,j,k) = 1/power(dy,2);
            at(i,j,k) = 1/power(dz,2); ab(i,j,k) = 1/power(dz,2);
            ap(i,j,k) = -ae(i,j,k)-aw(i,j,k)-an(i,j,k)-as(i,j,k)-at(i,j,k)-ab(i,j,k);
        end
    end
end

%calculate rhs term
for k=2:Nz+1
    for j=2:Ny+1
        for i=2:Nx+1
            rhs(i,j,k) = 50000*exp(-50*(power(1-x(i),2) + power(z(k),2)))*(100*(power(1-x(i),2) + power(z(k),2))-2);
        end
    end
end

%initialize phi value
for k=1:Nz+2
    for j=1:Ny+2
        for i=Nx+2
            phi(i,j,k) = 0;
            phi1(i,j,k) = 0;
        end
    end
end

L2norm = 1;
Tolerance = 0.00001;
count = 0;
%main loop
while (L2norm > Tolerance) && (count < 10000)
    %calculation of phi
    for k=2:Nz+1
        for j=2:Ny+1
            for i=2:Nx+1
                phi(i,j,k) = (rhs(i,j,k)-(ae(i,j,k)*phi(i+1,j,k)+aw(i,j,k)*phi(i-1,j,k)+an(i,j,k)*phi(i,j+1,k) ...
                    +as(i,j,k)*phi(i,j-1,k)+at(i,j,k)*phi(i,j,k+1)+ab(i,j,k)*phi(i,j,k-1)))/ap(i,j,k);
            end
        end
    end
    
    %update BC's
    for k=1:Nz+2
        for j=1:Ny+2
            for i=1:Nx+2
                phi(1,j,k) = 2*500*exp(-50*(1+power(z(k),2))) - phi(2,j,k);
                phi(Nx+2,j,k) = 2*(100*(1-z(k)) + 500*exp(-50*power(z(k),2))) - phi(Nx+1,j,k);
                phi(i,j,1) = 2*(100*x(i) + 500*exp(-50*power(1-x(i),2))) - phi(i,j,2);
                phi(i,j,Nz+2) = 2*500*exp(-50*(power(1-x(i),2)+1)) - phi(i,j,Nz+1);
                %periodic BC
                phi(i,1,k) = phi(i,Ny+1,k);
                phi(i,Ny+2,k) = phi(i,2,k);
            end
        end
    end
    
    %calculate L2norm
    sum = 0;
    for k=1:Nz+2
        for j=1:Ny+2
            for i=1:Nx+2
                sum = sum + power(phi(i,j,k)-phi1(i,j,k),2);
            end
        end
    end
    
    L2norm = sqrt(sum);
    count = count +1;
    
    %storing old values of phi
    for k=1:Nz+2
        for j=1:Ny+2
            for i=1:Nx+2
                phi1(i,j,k) = phi(i,j,k);
            end
        end
    end
    
    disp(count)
    disp(L2norm)
end

%Extracting values at y=Ny/2
for k=1:Nz+2
    for i=1:Nx+2
        u(k,i) = phi(i,Ny/2,k);
    end
end

%Analytical solution
for k=1:Nz+2
    for i=1:Nx+2
        v(k,i) = 500*exp(-50*(power(1-x(i),2)+power(z(k),2))) + 100*x(i)*(1-z(k));
    end
end

figure(1)
contourf(x,z,u,15)
colorbar
xlabel("X")
ylabel("Z")
title("GS for Nx=Ny=Nz=64")
figure(2)
contourf(x,z,v,15)
colorbar
xlabel("X")
ylabel("Z")
title("Analytical")
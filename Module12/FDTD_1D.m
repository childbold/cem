function [] = FDTD_1D(time_index)

% close all;

clc;

c = 3e8;          %Speed of light

dx = 1e-3;        %Spatial sampling

dt = dx/c;        %Temporal sampling

L = 1;            %Distance between plates in (m)

eps = 8.854e-12;

mu = 4*pi*1e-7;

Nx = round(L/dx) + 1;

Ez = zeros(Nx,1);

Hy = zeros(Nx,1);

Jz = zeros(Nx,1);

Nt = 2*time_index-1;

t = linspace(0,(Nt-1)*dt/2,Nt);

xe = linspace(0,1,Nx);

xh = linspace(dx/2,1-dx/2,Nx-1);

for n=1:length(t)
    
    if mod(n,2)==1
        
        Jz((Nx+1)/2) = exp(-((t(n)-2e-10)/5e-11)^2);
        
        for i=1:Nx-1
            
            Hy(i) = (dt/(mu*dx))*(Ez(i+1)-Ez(i))+Hy(i);
            
        end
        
    else
        for i=2:Nx-1
            
            Ez(i) = (dt/(eps*dx))*(Hy(i)-Hy(i-1))-(dt/eps)*Jz(i)+Ez(i);
            
        end
        
    end
    
end

figure(1)

    
plot(xe*1000,Ez,'r','LineWidth',2,'DisplayName','E_z');

hold on;

plot(xh*1000,377*Hy(1:Nx-1),'b--','LineWidth',2,'DisplayName','\eta_0H_y');

xlabel('x(mm)','FontSize',24)

ylabel('Field','FontSize',24)

legend('show')

grid on;

title(['Time step =',' ',num2str(time_index),',',' ','Time =',' ',num2str(t(Nt)/1e-9),' ',...
    'ns'],'FontSize',24);


end
        
        
        
        
    




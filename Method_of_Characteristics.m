%% Determination of Throat Section Area
clear; clc; close all;

gama = 1.4;
M = 2.4;
De = 2.3;   %Exit Diameter (m)
Ae = pi*De^2/4; %Exit Area (m^2)

syms At Dt
eq = (Ae / At)^2 == 1/(M^2) * ( 2/(gama+1) * ( 1 + (gama-1)/2 * M^2 ) )...
    ^((gama+1)/(gama-1));

At = double(vpasolve(eq,At));
At = At(2); %Throat Area (m^2)

eq2 = Dt^2*pi/4 == At;
Dt = double(vpasolve(eq2,Dt));
Dt = Dt(2); %Diameter of Throat (m)
save Dt

%% Divergent Section
clear; clc;

gama = 1.4; %Air Property
ncl = 10; %Number of Characteristics Lines

x(1) = 0;
y(1) = 0.4811;      %Change This Value to Achieve Exit Diameter
wx(1) = x(1);
wy(1) = y(1);

%Temporary Variables
l = ncl;
j = 1;

%Exit Conditions
M_exit = 2.4;

PM_exit = sqrt((gama + 1) / (gama-1) ) * atan (sqrt((gama-1)/(gama+1) * ...
    (M_exit^2-1) )) - atan(sqrt(M_exit^2-1));       %Radian
PM_exit = rad2deg(PM_exit);                 %Degree
teta_max = PM_exit / 2;                     %Degree


nnode = (2+(ncl+1)) * ncl/2;    %Number of nodes
node = (0:nnode)';

wall(ncl+1,1) = nnode;


for i = 1:(ncl-1)
    wall(ncl-i+1,1) = wall(ncl-i+2,1) - (i+1);   %Wall Points
end

teta_int = teta_max - floor(teta_max);  %Initial Turn Angle
inc = (teta_max - teta_int) / (ncl-1);  %Increment of Turn Angle

%Main Calculations
for i = 1 :nnode
    if i < wall(j+1)
        if wall(j) == 0 %Interactions With Point A
            teta(1) = teta_int + (i-1) * inc;
            PM_ang(1) = teta(1);
            Km(1) = teta(1) + PM_ang(1);
            Kp(1) = -Km(1);
            
            Mach(1,1) = Machfinder(PM_ang(1,1));
            mu(1,1) = asind(1/Mach(1,1));
            
            teta(i+1,1) = (Km(1) + Kp(i,1))/2;
            PM_ang(i+1,1) = (Km(1) - Kp(i,1))/2;
            Km(i+1,1) = teta(i+1,1) + PM_ang(i+1,1);
            Kp(i+1,1) = teta(i+1,1) - PM_ang(i+1,1);
            
            Mach(i+1,1) = Machfinder(PM_ang(i+1,1));
            mu(i+1,1) = asind(1/Mach(i+1,1));
            
            a_minus = (teta(1) + teta(i+1,1)) / 2 - ...
                (mu(1) + mu(i+1,1)) / 2;
            m_C_minus(i+1,1) = tand(a_minus);
            
            if i == 1   %For Point 1
                y(i+1,1) = 0;
                x(i+1,1) = -y(1) / (m_C_minus(i+1));
                m_C_plus(i+1,1) = -m_C_minus(i+1,1);
                centerlinex(i+1,1) = x(i+1,1);    
            else
                a_plus = (teta(i,1) + teta(i+1)) / 2 + (mu(i,1) + ...
                    mu(i+1)) / 2;
                m_C_plus(i+1,1) = tand(a_plus);
                x(i+1,1) = (y(1) - y(i) + m_C_plus(i+1) * x(i) - ...
                    m_C_minus(i+1) * x(1)) / (m_C_plus(i+1) - ...
                    m_C_minus(i+1));
                y(i+1,1) = y(i) + m_C_plus(i+1)* (x(i+1) - x(i));
            end     %Finish of A point Interations
            
        else %Grid Points Independent of A
            l = ncl - j + 1;
            Km(i+1,1) = Km(i-l);
            Kp(i+1) = - Km(ncl + 2 -l);
            teta(i+1,1) = (Km(i+1,1) + Kp(i+1,1)) / 2;
            PM_ang(i+1,1) = (Km(i+1,1) - Kp(i+1,1)) / 2;
            
            Mach(i+1,1) = Machfinder(PM_ang(i+1,1));
            mu(i+1,1) = asind(1/Mach(i+1,1));
            
            %Calculation of Slopes
            a_minus = (teta(i-l,1) + teta(i+1,1)) / 2 - (mu(i-l,1)...
                + mu(i+1,1)) / 2; 
            m_C_minus(i+1,1) = tand(a_minus);
            
            if teta(i+1) == 0  %Centerline
                y(i+1,1) = 0;
                x(i+1,1) = x(i-l) - y(i-l) / (m_C_minus(i+1));
                m_C_plus(i+1,1) = -m_C_minus(i+1,1);
                centerlinex(i+1,1) = x(i+1,1);
            else
                a_plus = (teta(i+1,1) + teta(i+1,1)) / 2 + (mu(i,1)...
                    + mu(i+1,1)) / 2;  
                m_C_plus(i+1,1) = tand(a_plus);
                x(i+1,1) = (y(i-l) - y(i) + m_C_plus(i+1) * x(i) - ...
                    m_C_minus(i+1) * x(i-l)) / (m_C_plus(i+1) - ...
                    m_C_minus(i+1));
                y(i+1,1) = y(i) + m_C_plus(i+1) * (x(i+1) - x(i));
            end
        end
        
    elseif  i == wall(j+1)    %Wall Points
        Km(i+1) = Km(1);
        Kp(i+1,1) = Kp(i,1);
        teta(i+1,1) = (Km(i+1,1) + Kp(i+1,1)) / 2;
        PM_ang(i+1,1) = (Km(i+1,1) - Kp(i+1,1)) / 2;
        
        Mach(i+1,1) = Machfinder(PM_ang(i+1,1));
        mu(i+1,1) = asind(1/Mach(i+1,1));
        a_wall = (teta(wall(j)+1) + teta(i)) / 2;
        a_plus = teta(i) + mu(i);
        
        m_C_wall(i+1,1) = tand(a_wall);
        m_C_plus(i+1,1) = tand(a_plus);
        x(i+1,1) = (y(wall(j)+1) - y(i) + m_C_plus(i+1,1) * x(i) - ...
            m_C_wall(i+1,1) * x(wall(j)+1)) / (m_C_plus(i+1,1) - m_C_wall(i+1)) ;
        y(i+1,1) = y(i) + m_C_plus(i+1,1) * (x(i+1) - x(i));
        
        wx(j+1,1) = x(i+1,1);
        wy(j+1,1) = y(i+1,1);
        
        j = j+1;
    end
end

m_C_minus(i+1) = 0;

%Display Outputs
Output.n(1) = "a ";
Output.n(2:nnode+1,1) = (node(2):nnode)';
Output.Km = Km;
Output.Kp = Kp;
Output.teta = teta;
Output.PM_ang = PM_ang;
Output.Mach = Mach;
Output.Machang = mu;
Output.slopem = m_C_minus;
Output.slopep = m_C_plus;
Output.walla = m_C_wall;

Output.x = x;
Output.y = y;

T = struct2table( Output );
T.Properties.VariableNames(1) = {'Point#'};
disp(T)

xi = linspace (min(x),max(x), 100);
yi_sphline = interp1(wx,wy,xi,'spline');
c(1:length(xi)) = 0;
centerlinex = nonzeros(centerlinex);

centerliney = zeros(length(centerlinex),1);

%%%%%%%% Plotting %%%%%%%
figure(1)
plot(xi,c,'--k'); hold on;
plot(centerlinex,centerliney,'.b'); hold on;
plot(wx,wy,'.b'); hold on;
plot(xi,yi_sphline,'k'); hold on;
for i = 1:length(centerlinex)
plot([wx(1) centerlinex(i)],[wy(1) centerliney(i)], '-r' ); hold on;
end

for i = 1:ncl
    plot( [centerlinex(i) wx(i+1)], [centerliney(i) wy(i+1)], '-b'  ); hold on
end

writetable(struct2table(Output), 'Project.xls')

%Function to Find Mach Number
f = @Machfinder;

function M = Machfinder(PM_ang)
gama = 1.4;
syms M
eq = - deg2rad(PM_ang) + sqrt((gama + 1) / (gama-1) ) * atan (sqrt((gama-1)/...
    (gama+1) * (M^2-1) )) - atan(sqrt(M^2-1)) == 0;
M = abs(double(vpasolve(eq,M)));

end
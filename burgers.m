clc
clear

%% QUEDA PER FER-
% - COMENTAR BÉ EL CODI
% - Penjar al GITHUB


%% 1. Comparing cases
clc
clear
% Parameters
Re = 40;
N = [20,100,20,20]; % Number of Fourier modes
C1 = 0.03;
Ck = [0.5,0,0.4523,0.05];
LES = [0,0,1,1];

% Calculations
for j = 1:length(N)
    [Ek{j},k{j},it,names{j},t{j}] = solve_case(N(j),C1,Ck(j),Re,LES(j));
    fprintf("1. Solved %s. Elapsed time is %.2f s\n", names{j},t{j})
end

% Plotting
figure(1)
for j = 1:length(N)
    loglog(k{j}, Ek{j},'o-','LineWidth',1,'MarkerSize',3);   
    xlabel('k (Wave number)')
    ylabel('Ek (Kinetic energy)')
    grid on
    hold on
end

% Slope m = -2
x_slope = 0:1:100;
y_slope = x_slope.^(-2);
names{end+1} = "Slope m=-2";
loglog(x_slope, y_slope, "--")
legend(names,"Location","best")

% Saving
saveas(figure(1),"Images/cases_comparison.png")
%% %%%%%%%%%%%%%%%%%%%% 2. Comparing Reynolds
clc
clear
% Parameters
Re = 10:20:80;
N = 100; 
C1 = 0.03;
Ck = 0;
LES = 0; % DNS

% Calculations
for j = 1:length(Re)
    [Ek{j},k{j},it,names{j},t{j}] = solve_case(N,C1,Ck,Re(j),LES);
    fprintf("2. Solved %s. Elapsed time is %.2f s\n", names{j},t{j})
end

% Plotting
figure(2)
for j = 1:length(Ek)
    loglog(k{j}, Ek{j},'o-','LineWidth',1,'MarkerSize',3);   
    xlabel('k (Wave number)')
    ylabel('Ek (Kinetic energy)')
    grid on
    hold on
end
% Slope m = -2
x_slope = 0:1:100;
y_slope = x_slope.^(-2);
names{end+1} = "Slope m=-2"
loglog(x_slope, y_slope, "--")
legend(names,"Location","best")

% Saving
saveas(figure(2),"reynolds_comparison_dns.png")

%% %%%%%%%%%%%%%%%%%%%%%%%% 3. Comparing N 
clc
clear
% Parameters
Re = 40;
N = 20:20:100; 
C1 = 0.03;
Ck = 0.4523;
LES = 1; % DNS

% Calculations
for j = 1:length(N)
    [Ek{j},k{j},it,names{j},t{j}] = solve_case(N(j),C1,Ck,Re,LES);
    fprintf("3. Solved %s. Elapsed time is %.2f s\n", names{j},t{j})
end

% Plotting
figure(3)
for j = 1:length(Ek)
    loglog(k{j}, Ek{j},'o-','LineWidth',1,'MarkerSize',3);   
    xlabel('k (Wave number)')
    ylabel('Ek (Kinetic energy)')
    grid on
    hold on
end

% Slope m = -2
x_slope = 0:1:100;
y_slope = x_slope.^(-2);
names{end+1} = "Slope m=-2";
loglog(x_slope, y_slope, "--")
legend(names,"Location","best")

% Saving
saveas(figure(3),"Images/n_comparison_les.png")

%% %%%%%%%%%%%%%%%%%%%%%%%% 4. Comparing Ck 
clc
clear
% Parameters
Re = 40;
N = 20; 
C1 = 0.03;
Ck = 0.05:0.2:1;
LES = 1; 

% Calculations
for j = 1:length(Ck)
    [Ek{j},k{j},it,names{j},t{j}] = solve_case(N,C1,Ck(j),Re,LES);
    fprintf("4. Solved %s. Elapsed time is %.2f s\n", names{j},t{j})
end

% Plotting
figure(4)
for j = 1:length(Ek)
    loglog(k{j}, Ek{j},'o-','LineWidth',1,'MarkerSize',3);   
    xlabel('k (Wave number)')
    ylabel('Ek (Kinetic energy)')
    grid on
    hold on
end

% Slope m = -2
x_slope = 0:1:100;
y_slope = x_slope.^(-2);
names{end+1} = "Slope m=-2"
loglog(x_slope, y_slope, "--")
legend(names,"Location","best")

% Saving
saveas(figure(4),"Images/ck_comparison.png")

%% %%%%%%%%%%%%%%%%%%%%%%%% 5. Comparing t c1
clc
clear
% Parameters
Re = 40;
N = 100; 
C1 = linspace(0.001, 0.3,4);
Ck = 0.4523;
LES = 0; % DNS

% Calculations
for j = 1:length(C1)
    [Ek{j},k{j},it,names{j},t(j)] = solve_case(N,C1(j),Ck,Re,LES);
    fprintf("5. Solved %s. Elapsed time is %.2f s\n", names{j},t(j))
end

%% %%%%%%%%%%%%%%
% Plotting 1
figure(5)
for j = 1:length(Ek)
    loglog(k{j}, Ek{j},'o-','LineWidth',1,'MarkerSize',3);   
    xlabel('k (Wave number)')
    ylabel('Ek (Kinetic energy)')
    grid on
    hold on
end

% Slope m = -2
x_slope = 0:1:100;
y_slope = x_slope.^(-2);
names{end+1} = "Slope m=-2";
loglog(x_slope, y_slope, "--")
legend(names,"Location","best")

% Saving
saveas(figure(5),"Images/t_cl_values_comparison_dns.png")

%% %%%%%%%%
% Plotting 2
figure(6)
plot(C1,t,'o-','LineWidth',1,'MarkerSize',3);   
xlabel('C1 parameter')
ylabel('Computation time [s]')
grid on
legend("DNS with N=100, Re=40.0, Ck=0.4523 ","Location","northeast")


% Saving
saveas(figure(6),"Images/t_c1_comparison_les.png")

%% %%%%%%%%%%%%%%%%%%%%%%%% 6. Comparing t N
clc
clear
% Parameters
Re = 40;
N = 10:25:100; 
C1 = 0.03;

Ck = 0.4523;
LES = 1; % DNS
for j = 1:length(N)
    [~,~,~,names{j},tL(j)] = solve_case(N(j),C1,Ck,Re,LES);
    fprintf("5. Solved %s. Elapsed time is %.2f s\n", names{j},tL(j))
end

Ck = 0;
LES = 0; % DNS
for j = 1:length(N)
    [~,~,~,names{j},tD(j)] = solve_case(N(j),C1,Ck,Re,LES);
    fprintf("7. Solved %s. Elapsed time is %.2f s\n", names{j},tD(j))
end

figure(7)
plot(N,tL,'o-','LineWidth',1,'MarkerSize',3);
hold on
plot(N,tD,'o-','LineWidth',1,'MarkerSize',3); 
xlabel('N')
ylabel('Computation time [s]')
grid on
legend("LES with C1=0.03, Re=40.0, Ck=0.4523 ", "DNS with C1=0.03, Re=40.0, Ck=0 ","Location","best")


% Saving
saveas(figure(7),"Images/t_n_comparison.png")

%% %%%%%%%%%%%%%%%%%%%%%%%% 7. Comparing t Re LES
clc
clear
% Parameters
Re = 10:25:100;
N = 40; 
C1 = 0.03;

Ck = 0.4523;
LES = 1; % DNS
for j = 1:length(Re)
    [~,~,~,names{j},tL(j)] = solve_case(N,C1,Ck,Re(j),LES);
    fprintf("5. Solved %s. Elapsed time is %.2f s\n", names{j},tL(j))
end

N = 100; 
Ck = 0;
LES = 0; % DNS
for j = 1:length(Re)
    [~,~,~,names{j},tD(j)] = solve_case(N,C1,Ck,Re(j),LES);
    fprintf("7. Solved %s. Elapsed time is %.2f s\n", names{j},tD(j))
end

figure(8)
plot(Re,tL,'o-','LineWidth',1,'MarkerSize',3);
hold on
plot(Re,tD,'o-','LineWidth',1,'MarkerSize',3); 
xlabel('Re')
ylabel('Computation time [s]')
grid on
legend("LES with C1=0.03, N = 20, Ck=0.4523 ", "DNS with C1=0.03, N = 100, Ck=0 ","Location","best")


% Saving
saveas(figure(8),"Images/t_re_comparison.png")

%% %%%%%%%%%%%%%%%%%%%%%%%% 8. Comparing t Re DNS
clc
clear
% Parameters
Re = 5:25:150;
N = 100; 
C1 = 0.03;
Ck = 0;
LES = 0; % DNS

% Calculations
for j = 1:length(Re)
    [Ek{j},k{j},it,names{j},t(j)] = solve_case(N,C1,Ck,Re(j),LES);
    fprintf("5. Solved %s. Elapsed time is %.2f s\n", names{j},t(j))
end

%% %%%%%%%%
% Plotting 2
figure(8)
plot(Re,t,'o-','LineWidth',1,'MarkerSize',3);   
xlabel('Re')
ylabel('Computation time [s]')
grid on
legend("DNS with C1=0.03, N=100.0, Ck=0 ","Location","best")


% Saving
saveas(figure(8),"Images/t_re_comparison_dns.png")
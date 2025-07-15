%plot model trends compared with S&S Table 1

SS_data = [3520	4670	5370	5600	5560	4326	5530	5370	5060	4500	5370	4830	4160; %Ts
    0.33	4.57	8.81	12.2	14.4	2.16	7.96	8.81	8.24	6.4	8.81	5.39	2.77; %
    96	92	87	81	74	96	89	87	87	86	87	92	94] %



%pressure amplitude, bar
PA = [1.1	1.15	1.2	1.25	1.3	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2]
%init radius, um
R0 = [4.5	4.5	4.5	4.5	4.5	2	3	4.5	6.5	10	4.5	4.5	4.5]
%frequency, kHz
omega = [26.5	26.5	26.5	26.5	26.5	26.5	26.5	26.5	26.5	26.5	26.5	38.3	58.9] 

y_lab = 'Temperature (K)'


figure('Position',[400 400 800 260]);
subplot(1,3,1)
box on
hold on

plot(PA(1:5),out(1:5,1)','o',MarkerEdgeColor="#0072BD",MarkerFaceColor="#0072BD")

plot(PA(1:5),SS_data(1,1:5),'^',MarkerEdgeColor="#D95319",MarkerFaceColor="#D95319")

hold off
xlim([1,1.4])
ylim([0,6000])
xlabel('Pressure amplitude (bar)')
ylabel(y_lab)

%% 
subplot(1,3,2)
box on
hold on

plot(R0(6:10),out(6:10,1)','o',MarkerEdgeColor="#0072BD",MarkerFaceColor="#0072BD")

plot(R0(6:10),SS_data(1,6:10),'^',MarkerEdgeColor="#D95319",MarkerFaceColor="#D95319")

hold off
xlim([0,12])
ylim([0,6000])
xlabel('Initial radius ({\mu}m)')
ylabel(y_lab)


%%
subplot(1,3,3)
box on
hold on

plot(omega(11:end),out(11:end,1)','o',MarkerEdgeColor="#0072BD",MarkerFaceColor="#0072BD")

plot(omega(11:end),SS_data(1,11:end),'^',MarkerEdgeColor="#D95319",MarkerFaceColor="#D95319")

hold off
xlim([0,100])
ylim([0,6000])
xlabel('Frequency (kHz)')
ylabel(y_lab)

f = gcf;
exportgraphics(f,'fig_S1_temp.png','Resolution',300)

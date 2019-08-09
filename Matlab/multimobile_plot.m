% Plot the theoretical curve for the simulation with multiple antennas on
% mobiles

hold on
% lambda=num_mobile/(len_all^2*pi);
antenna_r_alpha=antenna_r.^-alpha;
pk=zeros(1,length(l_actual));
for i=1:length(l_actual)
pk(i)=l_actual(i)*sum(antenna_r_alpha);
end
eta=zeros(1,length(l_actual));
for l=1:length(l_actual)
eta(l)=M*log2(1+pk(l).*(alpha/(2*pi^2*lambda*M)*sin(2*pi/alpha)).^(alpha/2)*(K*l_actual(l)).^(alpha/2-1));
end
semilogx(l_actual,eta,'r*--')
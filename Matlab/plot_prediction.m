% Plot the theoretical results of the celluar model

% lambda=num_mobile/(len_all^2*pi);
antenna_r_alpha=antenna_r.^-alpha;
pk=zeros(1,length(antenna_r_alpha));
for i=1:length(antenna_r_alpha)
pk(i)=sum(antenna_r_alpha(1:i));
end
eta=zeros(1,length(k_actual));
for k=1:length(k_actual)
eta(k)=log2(pk(k_actual(k)).*(alpha/(2*pi^2*lambda)*sin(2*pi/alpha)).^(alpha/2)*(k_actual(k)).^(alpha/2-1));
end
semilogx(k_actual,eta,'g*--')


% A simulation comparing cellular and cell free model

% parameters
% k = number of receiving antennas
% len_all = side length of the big square "world"
% len_ctr = side length of the center square where antennas are deployed
% lambda = density of interfering mobile
% num_mobile = number of interfering mobile
% alpha = path loss exponent
rng('default');
rng(1);
len_all=10000;
len_ctr=100;
lambda=5*10^(-6);
num_mobile=len_all^2*pi*lambda;
alpha=4;    
% represent antenna locations as a complex number
k_temp=10000;
% assume the targete mobile is at 0+0i
antenna_r = 1+(len_ctr-1)*sqrt(rand(k_temp,1));
antenna_pos = antenna_r.*exp(2*pi*rand(k_temp,1)*1i);
reps=10;

% Create a semilog figure
% figure
% semilogx(0,0)
% hold on

k_actual=round(logspace(0,3,20));
a_mc_buffer=zeros(length(k_actual),reps);
a_bs_buffer=zeros(length(k_actual),reps);
for k=1:length(k_actual)
    for j=1:reps
        % calculate channel vector for target
        h1=(randn(k_actual(k),1)+1i*randn(k_actual(k),1)).*antenna_r(1:k_actual(k)).^(-alpha/2);

        % construct a channel coefficient matrix 
        H_cof=zeros(k_actual(k),num_mobile+1);
        H_cof(:,1)=h1;

        itf_pos = (1+(len_all-1)*sqrt(rand(num_mobile,1))).*exp(2*pi*rand(num_mobile,1)*1i);
        itf_r = zeros(k_actual(k),num_mobile);
        for i=1:num_mobile
            diff_itf=antenna_pos(1:k_actual(k))-itf_pos(i);
            itf_r(:,i)=abs(diff_itf);
        end

        for i=1:num_mobile
            h_i=(randn(k_actual(k),1)+1i*randn(k_actual(k),1)).*itf_r(:,i).^(-alpha/2);
            H_cof(:,i+1)=h_i;
        end

        % Do not include h1 for MMSE's SNR calculation
        R=H_cof(:,2:end)*H_cof(:,2:end)';
        a_mc_buffer(k,j)=real(h1'*(R\h1));
%         hold on
    end
end
% semilogx(k_actual,a_mc_buffer,'b+')
% xlabel('Number of antennas, K')
% ylabel('Spectral Efficiency')
% title('Cell-free(Blue) vs Single-tower(Red) Simulation')
% hold on
% grid on

for k=1:length(k_actual)
    for j=1:reps
        % calculate channel vector for target
        d=sum(antenna_r(1:k_actual(k)).^(-alpha))/k_actual(k);
        h1=(randn(k_actual(k),1)+1i*randn(k_actual(k),1)).*d.^(1/2);

        % construct a channel coefficient matrix 
        H_cof=zeros(k_actual(k),num_mobile+1);
        H_cof(:,1)=h1;

%         itf_pos = (1+(len_all-1)*rand(num_mobile,1)).*exp(2*pi*rand(num_mobile,1)*1i);
        itf_r = 1+(len_all-1)*sqrt(rand(num_mobile,1));

        for i=1:num_mobile
            h_i=(randn(k_actual(k),1)+1i*randn(k_actual(k),1)).*itf_r(i).^(-alpha/2);
            H_cof(:,i+1)=h_i;
        end

        % Do not include h1 for MMSE's SNR calculation
        R=H_cof(:,2:end)*H_cof(:,2:end)';
        a_bs_buffer(k,j)=real(h1'*(R\h1));
%         semilogx(k_actual(k),a,'r.')
%         hold on
    end
end
% semilogx(k_actual,a_bs_buffer,'r.')
% 
% hold off
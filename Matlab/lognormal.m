% A simulation with log-normal fading

% parameters
% k = number of receiving antennas
% len_all = radius of the big square "world"
% len_ctr = radius of the center square where antennas are deployed
% lambda = density of interfering mobile
% num_mobile = number of interfering mobile
% alpha = path loss exponent
% M = the number of antennas per mobile
% L = the number of antennas per base station
% sigma = variance for lognormal variable
rng('default');
rng(1);
len_all=100000;
len_ctr=100;
lambda=5*10^(-6);
num_mobile=len_all^2*pi*lambda;
alpha=4;    
% represent antenna locations as a complex number
K=2;
L_temp=100;
M=1;
sigma=8;

% assume the targete mobile is at 0+0i

antenna_r = 1+(len_ctr-1)*sqrt(rand(K,1));
bs_pos = antenna_r.*exp(2*pi*rand(K,1)*1i);
reps=100;

% Create a semilog figure
figure
l_actual=round(logspace(0,3,20));
a_mc_buffer=zeros(length(l_actual),reps);

% Make a fixed set of lognormal random numbers for the channel of target
mySample=10.^(randn(K,1)*sigma/10);
for l=1:length(l_actual)
    for j=1:reps
        % calculate channel vector for target
        h1=(randn(l_actual(l)*K,M)./sqrt(2)+1i*randn(l_actual(l)*K,M)./sqrt(2));
        for k=1:K
            h1(1+l_actual(l)*(k-1):l_actual(l)*k,:)=h1(1+l_actual(l)*(k-1):l_actual(l)*k,:).*antenna_r(k).^(-alpha/2).*mySample(k).^(1/2);
        end  
        
        % construct a channel coefficient matrix 
        H_cof=randn(K*l_actual(l), M*num_mobile)./sqrt(2)+1i*randn(K*l_actual(l), M*num_mobile)./sqrt(2);
%         H_cof(:,1)=h1;
        for i=1:num_mobile
            for k=1:K
                H_cof(1+l_actual(l)*(k-1):k*l_actual(l),i)=H_cof(1+l_actual(l)*(k-1):k*l_actual(l),i).*(10.^(randn(l_actual(l),1)*sigma/10)).^(1/2);
            end
        end

        itf_pos = (1+(len_all-1)*sqrt(rand(num_mobile,1))).*exp(2*pi*rand(num_mobile,1)*1i);
        % Each column 
        itf_r = zeros(K,num_mobile);
        diff_itf = zeros(num_mobile,1);
        for i=1:K
            diff_itf=bs_pos(i)-itf_pos;
            itf_r(i,:)=abs(diff_itf);
        end
        
        for k=1:K
            for m=1:num_mobile
                % this doesn't look like it's multiplying
                H_cof(1+l_actual(l)*(k-1):l_actual(l)*k,1+M*(m-1):M*m)=H_cof(1+l_actual(l)*(k-1):l_actual(l)*k,1+M*(m-1):M*m).*itf_r(k,m).^(-alpha/2);
            end
        end

        % Do not include h1 for MMSE's SNR calculation
        R=H_cof(:,2:end)*H_cof(:,2:end)';
        a_mc_buffer(l,j)=log2(det(eye(M)+real(h1'*(R\h1))));
    end
end
semilogx(l_actual,a_mc_buffer,'b+')
xlabel('Number of antennas per base station, L')
ylabel('Spectral Efficiency')
title('Multiple Antennas on Mobiles')
grid on

hold off
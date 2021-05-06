clc
clear
close all
%% head
% @param choice_pilots: Pivots sets. Values are [A,B] 
% @param choice_channel: 1/4/8-tap multipath channel. Values are [A,B,C] 
% @param subc: select subcarriers, values are in range(1,64) 
%              nclusive. [10,25] for question prompt
% @param snr: snr. [5,10]db for question prompt

choice_pilots = 'B'; % select pilots set
choice_channel = 'A'; % select channel
subc = 10; % subcarriers
snr = 10; % select snr
N = 5e3; % num of symbols


%% body
% Init
M = 4; % Modulation order
K = 64; % Number of sub-carriers
b_in = randi([0 M-1], N, K)' ;
%
h1 = 1;
h2 = [-0.7298 - 0.2884i, -0.1235 + 0.0411i, -0.0108 - 0.0381i, -0.0654 - 0.0165i];
h3 = [-1.4896 + 0.7091i, 0.4849 + 0.1602i, 0.0339 - 0.0344i, -0.2158 + 0.0961i,...
- 0.0675 + 0.0593i, 0.1587 - 0.0241i, -0.0049 - 0.0335i, -0.0232 + 0.0434i];

sub_car1 = 1;
sub_car2 = [8, 24, 40, 56];
sub_car3 = [4, 12, 20, 28, 36, 44, 52, 60];

% 4 QAM
X = qammod(b_in,M) ./ sqrt(2) ; 

% add pilots
if choice_pilots == 'A'
    start = 8;
    while start<=K
        X(start,:) = 5+5i;
        start = start + 16;
    end
end
if choice_pilots == 'B'
    start = 4;
    while start<=K
        X(start,:) = 1+1i;
        start = start + 8;
    end
end
% select multipath
if choice_channel == 'A'
    h = h1;
    num_multipath = 1;
    car = sub_car1;
end
if choice_channel == 'B'
    h = h2;
    num_multipath = 4;
    car = sub_car2;
end
if choice_channel == 'C'
    h = h3;
    num_multipath = 8;
    car = sub_car3;
end

%FFT
x = ifft(X,K);
% P2S + cyclic prefix
x_in = zeros(K+16,N);
for i = 1 : N
    prefix = x(end-15:end,i);
    data = x(:,i);
    x_in(:,i) = vertcat(prefix,data);
end

%% Channel
for i = 1 : N
    cur_symbol = conv(h,x_in(:,i));
    pwr_avg = sum(cur_symbol.*conj(cur_symbol))/length(cur_symbol);
    pwr_avg_db = pow2db(pwr_avg);
    % adjust snr
    y_out(:,i) = awgn(cur_symbol(1:80),snr,pwr_avg_db,'dB');
end


%% Receiver

y = y_out(end-63:end,:);
Y = fft(y);

% hyper param for choosing sub carrier



X_hat = zeros(64,1);
for n = 1 : N
    j = 1;
    counter = 0;
    for i = 1 : K
        counter = counter + 1;

        h = car(j);
        H = conj(X(h,n)) * Y(h,n) / (abs(X(h,n)))^2;
        X_hat(i,n) = conj(H)*Y(i,n)/(abs(H))^2;

        if counter == 64/num_multipath
            counter = 0;
            j = j+1;
        end
    end
end

constDiagram = comm.ConstellationDiagram();
constDiagram(X_hat(subc,:)')
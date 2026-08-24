clc,clear
a = '++---+-----+++-+++++-++-++-+-+-++--+++-+--+---+++-++--+-++-+---+++------++';  % 序列 A
b = '-+++-+-+--++-++---+++++-++++-+-++++--++--++-+--+--+++-+-+------+-++++-++++';  % 序列 B

% ------------------------------------------
% 以下程式碼會自動計算週期性自相關 (Periodic ACF)
% ------------------------------------------

% 1. 將符號轉換為數值 (+ -> 1, - -> -1)
s1 = 2 * double(a == '+') - 1;
s2 = 2 * double(b == '+') - 1;

% 檢查長度是否一致 (通常互補對長度需相同)
L = length(s1);
if length(s2) ~= L
    warning('注意：兩個序列的長度不一致！計算可能無意義。');
end

% 2. 計算週期性自相關函數 (Periodic Autocorrelation)
%    定義: rho(u) = sum( s(k) * s(k+u) )，其中 k+u 為循環位移
pacf1 = zeros(1, L);
pacf2 = zeros(1, L);

for u = 0 : L-1
    pacf1(u+1) = sum(s1 .* circshift(s1, -u));
    pacf2(u+1) = sum(s2 .* circshift(s2, -u));
end

% 3. 計算兩者總和 (用於檢查是否為 PCP)
pacf_sum = pacf1 + pacf2;

% 4. 輸出數值結果到 Command Window
fprintf('\n================ 週期性相關結果 ================\n');
fprintf('序列長度 L = %d\n', L);
fprintf('----------------------------------------------\n');

% 顯示序列 a 結果
fprintf('1. 序列 a (%s) 的週期自相關:\n', a);
fprintf('Lag: '); fprintf('%3d ', 0:L-1); fprintf('\n');
fprintf('Val: '); fprintf('%3d ', pacf1); 
fprintf('\n\n');

% 顯示序列 b 結果
fprintf('2. 序列 b (%s) 的週期自相關:\n', b);
fprintf('Lag: '); fprintf('%3d ', 0:L-1); fprintf('\n');
fprintf('Val: '); fprintf('%3d ', pacf2); 
fprintf('\n\n');

% 顯示總和結果
fprintf('3. 兩序列週期自相關函數之和 (Sum):\n');
fprintf('Lag: '); fprintf('%3d ', 0:L-1); fprintf('\n');
fprintf('Val: '); fprintf('%3d ', pacf_sum); 
fprintf('\n');

% 5. 自動判斷是否為週期性互補對 (PCP)
%    檢查除了 Lag 0 以外的所有值是否為 0
sidelobes = pacf_sum(2:end); % 取出 Lag 1 到 L-1

if all(sidelobes == 0)
    fprintf('\n[判定] YES! 這是一組週期性互補對 (PCP)。\n');
    fprintf('       (Lag 0 的值為 %d，其餘皆為 0)\n', pacf_sum(1));
else
    fprintf('\n[判定] NO。這不是 PCP，旁瓣不全為 0。\n');
end
fprintf('==============================================\n');

% ------------------------------------------
% 以下程式碼計算 PAPR (Peak-to-Average Power Ratio)
% ------------------------------------------
fprintf('\n正在計算 PAPR...\n');

% 設定過取樣倍率 (Oversampling Factor)
% 通常取 4 或 8 即可獲得足夠精確的連續時間 PAPR 近似值
ovs = 8; 
N_fft = L * ovs; % FFT 點數

% 1. 計算過取樣 FFT (將序列補零至 N_fft 長度)
F1 = fft(s1, N_fft);
F2 = fft(s2, N_fft);

% 2. 計算功率頻譜 (Power Spectrum)
P1 = abs(F1).^2;
P2 = abs(F2).^2;

% 3. 計算 PAPR (線性值)
% 定義: PAPR = max(Power) / Average(Power)
% 對於 {-1, 1} 二元序列，理論平均功率 (Average Power) 等於序列長度 L
% 這裡我們使用 mean(P) 來計算平均值，結果應約等於 L
papr1_linear = max(P1) / mean(P1);
papr2_linear = max(P2) / mean(P2);

% 4. 轉換為分貝 (dB)
papr1_dB = 10 * log10(papr1_linear);
papr2_dB = 10 * log10(papr2_linear);

% 5. 輸出結果
fprintf('\n================ PAPR 計算結果 ================\n');
fprintf('過取樣倍率 (Oversampling): %d\n', ovs);
fprintf('----------------------------------------------\n');
fprintf('1. 序列 a 的 PAPR:\n');
fprintf('   Linear: %8.4f\n', papr1_linear);
fprintf('   dB    : %8.4f dB\n', papr1_dB);
fprintf('\n');
fprintf('2. 序列 b 的 PAPR:\n');
fprintf('   Linear: %8.4f\n', papr2_linear);
fprintf('   dB    : %8.4f dB\n', papr2_dB);
fprintf('==============================================\n');
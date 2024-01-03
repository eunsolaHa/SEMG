clc; close all; clear;

%% Surface electromyography signal analysis
  % �ð� ���� : �������������(RMS)
  % �ð�-���ļ� ���� : ������ļ�(MNF), �߾����ļ�(MDF)

% folder�� �ִ� '.txt' ���� �ҷ�����
folder = 'C:\Users\jang\Desktop\������ ������\data';

files = dir(fullfile(folder, '*.txt'));

for file_1 = 1:length(files)
    data = load(fullfile(folder, files(file_1).name));
    xin = data(:,1);
    
    %sEMG �ð�ȭ
    figure;
    plot(xin)
    xlabel('Time(ms)');
    ylabel('Voltage(mV)');
    title('sEMG Raw Data');
    
    [filepath, filename, ~ ] = fileparts(files(file_1).name);
    
    % ���� ����
    save_folder1 = fullfile(folder, 'RawEMG');
    if ~exist(save_folder1, 'dir')
        mkdir(save_folder1);
    end
    
    % ���� ���� ����
    %save_filename_xlsx = fullfile(save_folder1, [filename, '.xlsx']);
    %T1 = table(xin, 'VariableNames', {'Raw_EMG_signal'});
    %writetable(T1, save_filename_xlsx, 'Sheet', 1);
    
    % PNG ���� ����
    %save_filename_png = fullfile(save_folder1, [filename, '.png']);
    %saveas(gcf, save_filename_png);
    
    %4�� butterwoth bandpass filter
    fs = 1000;
    fpass = [20, 450];
    
    [b, a] = butter(4, [fpass(1)/ (fs/2) fpass(2)/(fs/2)], 'bandpass');
    [H W] = freqz(b,a,1000);
    
    BBF_xin = filter(b,a, xin);
    
    figure;
    p = plot(BBF_xin, 'b');
    set(p, 'LineWidth', 0.1);
    xlabel('Time(ms)');
    ylabel('Voltage(mV)');
    title('Filtered sEMG Data'); 
    
   % ���� ����
    %save_folder2 = fullfile(folder, 'FilteredEMG');
    %if ~exist(save_folder2, 'dir')
    %    mkdir(save_folder2);
    %end
    
    % ���� ���� ����
    %save_filename_xlsx = fullfile(save_folder2, [filename, '.xlsx']);
    %T2 = table(xin, 'VariableNames', {'Filtered_EMG_signal'});
    %writetable(T2, save_filename_xlsx, 'Sheet', 1);
    
    % PNG ���� ����
    %save_filename_png = fullfile(save_folder2, [filename, '_filtered.png']);
    %saveas(gcf, save_filename_png);
    
    
    % 4th order bandpass characteristic output
    %digital filter frequency W -> real frequency vector w
    
    w = W * (fs/2) / pi;
    plot(w, abs(H)); axis([0 ,max(w) 0 1.2]);
    xlabel('Frequency(Hz)');
    ylabel('Voltage(mV)');
    title('4th order Butterworth bandpass filter');
    
    % ���� ����
    %save_folder3 = fullfile(folder, '4bandpassfilter');
    %if ~exist(save_folder3, 'dir')
    %    mkdir(save_folder3);
    %end
    
    % ���� ���� ����
    %save_filename_xlsx = fullfile(save_folder3, [filename, '_filtered.xlsx']);
    %T3 = table(xin, 'VariableNames', {'Filtered_EMG_signal'});
    %writetable(T3, save_filename_xlsx, 'Sheet', 1);
    
    % PNG ���� ����
    %save_filename_png = fullfile(save_folder3, [filename, '_filtered.png']);
    %saveas(gcf, save_filename_png);
    
    % rolling window
    
    l = length(BBF_xin); %l = window Number
    window = 500; % window = window length
    r = floor(l/window); % r = window count
    z = zeros(window, r);
    
    for ii = 1:r
        start = int32(500*(ii-1)+1);
        finish = int32(500*ii);
        z(:, ii) = BBF_xin(start:finish,1);
    end
    
    %RMS
    rms_c = zeros(r,1);
    for iii = 1:r
        rms_c(iii,1) = rms(z(:,iii));
        %rms_c(iii, 1) = sqrt(mean(z(:, iii).^2));
    end
    
    figure;
    plot(rms_c)
    xlabel('Time');
    ylabel('Amplitude');
    title('RMS(�������������)');
    
    % ���� ����
    %save_folder4 = fullfile(folder, 'RMS');
    %if ~exist(save_folder4, 'dir')
        %mkdir(save_folder4);
    %end
    
    % ���� ���� ����
    %save_filename_xlsx = fullfile(save_folder4, [filename, '_RMS.xlsx']);
    %T4 = table(xin, 'VariableNames', {'RMS'});
    %writetable(T4, save_filename_xlsx, 'Sheet', 1);
    
    % PNG ���� ����
    %save_filename_png = fullfile(save_folder4, [filename, '_RMS.png']);
    %saveas(gcf, save_filename_png);
    
    %RMS ���� Ȯ�� ���� �ð�ȭ
    figure;
    histogram(rms_c);
    xlabel('Amplitude');
    ylabel('Frequency');
    title('���� ����');
    
    % ���� ����
    %save_folder5 = fullfile(folder, 'RMS_histogram');
    %if ~exist(save_folder5, 'dir')
    %    mkdir(save_folder5);
    %end
    
    % ���� ���� ����
    %save_filename_xlsx = fullfile(save_folder5, [filename, '_RMS_hist.xlsx']);
    %T5 = table(xin, 'VariableNames', {'RMS_Amplitude'});
    %writetable(T5, save_filename_xlsx, 'Sheet', 1);
    
    % PNG ���� ����
    %save_filename_png = fullfile(save_folder5, [filename, '_RMS_hist.png']);
    %saveas(gcf, save_filename_png);
    
    emgSignal = BBF_xin
    fs = 1000;
    windowLength = 500;
    gaussian = gausswin(windowLength); % gaussian window
    Noverlap = 0
    [s, f, t] = spectrogram(BBF_xin, gaussian, Noverlap, windowLength, fs);
    
    figure;
    imagesc(t, f, abs(s));
    axis xy;
    xlabel('Time(s)');
    ylabel('Frequency(Hz)');
    title('STFT(short-time fourier transform) ');
    %xlim([0 10])
    %ylim([0 max(f)]);
    
    %Mean frequency
    mnf_f = meanfreq(abs(s), f);
    %mnf_f(1) = mean(mnf_f);
    
    figure;
    plot(t, mnf_f) % �ð��� ���� ��� ���ļ� 
    plot(mnf_f) % �ð� ���� ���� �� �ð� �����̽��� ��� ���ļ�
    xlabel('time(s)');
    ylabel('Frequency (Hz)');
    title('������ļ�(Mean Frequency)');
    
    % ���� ����
    %save_folder6 = fullfile(folder, 'MNF');
    %if ~exist(save_folder6, 'dir')
     %   mkdir(save_folder6);
    %end
    
    % ���� ���� ����
    %save_filename_xlsx = fullfile(save_folder6, [filename, '_MNF.xlsx']);
    %T6 = table(xin, 'VariableNames', {'MeanFrequency'});
    %writetable(T6, save_filename_xlsx, 'Sheet', 1);
    
    % PNG ���� ����
    %save_filename_png = fullfile(save_folder6, [filename, '_MNF.png']);
    %saveas(gcf, save_filename_png);
    
    %MDF
    mdf_f = medfreq(abs(s), f); % �ð��� ���� �߾����ļ��� ��ȭ
    %mean_mdf = mean(mdf_f)
    
    figure;
    plot(t, mdf_f);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('�߾����ļ�(MDF)');
    
    % ���� ����
    %save_folder7 = fullfile(folder, 'MDF');
    
    %if ~exist(save_folder7, 'dir')
    %    mkdir(save_folder7);
    %end
    
    % ���� ���� ����
    %save_filename_xlsx = fullfile(save_folder7, [filename, '_mdf.xlsx']);
    %T7 = table(xin, 'VariableNames', {'Median_Frequency'});
    %writetable(T7, save_filename_xlsx, 'Sheet', 1);
    
    % PNG ���� ����
    %save_filename_png = fullfile(save_folder7, [filename, '_mdf.png']);
    %saveas(gcf, save_filename_png);
    
    
    
    %MDF �߼��� �߰�
    mdf_f = medfreq(abs(s), f); % �ð��� ���� �߾����ļ��� ��ȭ
    mean_mdf = mean(mdf_f)
    
    % �߼��� ��� (���׽� ������ �����Ͽ� �ʿ信 ���� ����)
    p = polyfit(t, mdf_f, 1);
    trendline = polyval(p, t);
    
    figure;
    plot(t, mdf_f, 'b', t, trendline, 'r--');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    legend('�߾����ļ�', '�߾����ļ� ����');
    title('� �ð��� ���� �߾����ļ� ��ȭ');
    xlim([0 54])
 
    
    % �߼����� ���� ���׽��� �׷����� ǥ��
    equation_text = sprintf('Y: %.4f * t + %.4f', p(1), p(2));
    text(0.5, max(mdf_f), equation_text, 'Color', 'r');
    text(mean(t), 1.05 * max(mdf_f), equation_text, 'Color', 'r', 'HorizontalAlignment', 'center');
    
    % ���� ����
    %save_folder8 = fullfile(folder, 'MDF_trend');
    %if ~exist(save_folder8, 'dir')
    %    mkdir(save_folder8);
    %end
    
    % ���� ���� ����
    %save_filename_xlsx = fullfile(save_folder8, [filename, '_MDF.xlsx']);
    %T8 = table(xin, 'VariableNames', {'MDF'});
    %writetable(T8, save_filename_xlsx, 'Sheet', 1);
    
    % PNG ���� ����
    %save_filename_png = fullfile(save_folder8, [filename, '_MDF.png']);
    %saveas(gcf, save_filename_png);
    
   
end

    

close all; clear all; clc

data_1 = load('Min_Fitness.txt');
data_2 = load('Ave_Fitness.txt');
data_3 = load('Max_Fitness.txt');


num_gen = 300;
num_sr = 30;
sum_data_1 = zeros(1, num_gen);
sum_data_2 = zeros(1, num_gen);
sum_data_3 = zeros(1, num_gen);
err_data_1 = zeros(1, num_gen);
err_data_2 = zeros(1, num_gen);
err_data_3 = zeros(1, num_gen);
Generation = 1:1:num_gen;

%Sums each gen from each sr
for (sr=1:num_sr)
    for (gen=1:num_gen)
        sum_data_1(gen) = sum_data_1(gen) + data_1(sr,gen);
        sum_data_2(gen) = sum_data_2(gen) + data_2(sr,gen);
        sum_data_3(gen) = sum_data_3(gen) + data_3(sr,gen);
    end
end
%Calculates the average
Ave_data_1 = sum_data_1/num_sr;
Ave_data_2 = sum_data_2/num_sr;
Ave_data_3 = sum_data_3/num_sr;

%Calculates the error bars for each generation
for (gen=1:num_gen)
    stdv_ph_1 = zeros(1, num_sr);
    stdv_ph_2 = zeros(1, num_sr);
    stdv_ph_3 = zeros(1, num_sr);
    for (sr=1:num_sr)
        stdv_ph_1(sr) = data_1(sr,gen);
        stdv_ph_2(sr) = data_2(sr,gen);
        stdv_ph_3(sr) = data_3(sr,gen);
    end
    if (~mod(gen,20))
    err_data_1(gen) = std(stdv_ph_1)/sqrt(num_sr);
    err_data_2(gen) = std(stdv_ph_2)/sqrt(num_sr);
    err_data_3(gen) = std(stdv_ph_3)/sqrt(num_sr);
    elseif(gen == 1)
    err_data_1(gen) = std(stdv_ph_1)/sqrt(num_sr);
    err_data_2(gen) = std(stdv_ph_2)/sqrt(num_sr);
    err_data_3(gen) = std(stdv_ph_3)/sqrt(num_sr);
    else
    err_data_1(gen) = NaN;
    err_data_2(gen) = NaN;
    err_data_3(gen) = NaN;
    end
end

figure(1)
hold all;
errorbar(Generation, Ave_data_1, err_data_1, 'LineWidth', 1.5)
errorbar(Generation, Ave_data_2, err_data_2, 'LineWidth', 1.5)
errorbar(Generation, Ave_data_3, err_data_3, 'LineWidth', 1.5)
title('DeJong 1')
xlabel('Generation', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Fitness', 'FontSize', 20, 'FontWeight', 'bold')
set(gca,'FontSize',15)
axis([0 num_gen 0 100])
legend ('min', 'ave', 'max')
%print('R1','-dpng')
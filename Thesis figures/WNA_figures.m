%% Raster plot
dataDirectory = "2022.10.14";
filename = "M1_N1_blwgn";

P = getStructP(dataDirectory, filename,[nan nan],1);
plot_data(P(1), "stimulus", "raster", "gcfr");

%% Jitter histogram
c = parula(height(T_jitter_fidelity));
figure; colororder(c);
for i = 1:height(T_jitter_fidelity)
histogram(cell2mat(T_jitter_fidelity.jitter(i)), 'BinWidth', 0.1,...
'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'Normalization','probability',...
'DisplayName',replace(join([T_jitter_fidelity.date(i) T_jitter_fidelity.filename(i)]," "), "_", " "));
hold on;
end
set(gca, 'YLimitMethod', 'padded');
xlim([-0.25 4]);
box off; ylabel('Probability'); xlabel('Jitter (ms)');
mean_jitter = mean(cell2mat(T_jitter_fidelity.jitter));
std_jitter = std(cell2mat(T_jitter_fidelity.jitter));

text_str = sprintf("Mean \\pm STD = %0.2f \\pm %0.2f ms", mean_jitter, std_jitter);
text(2,0.4, text_str);

%% Fidelity histogram
c = parula(height(T_jitter_fidelity));
figure; colororder(c);
for i = 1:height(T_jitter_fidelity)
histogram(cell2mat(T_jitter_fidelity.fidelity(i)), 'BinWidth', 0.05,...
'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'Normalization','probability',...
'DisplayName',replace(join([T_jitter_fidelity.date(i) T_jitter_fidelity.filename(i)]," "), "_", " "));
hold on;
end
set(gca, 'YLimitMethod', 'padded');
xlim([-0.05 1]);
box off; ylabel('Probability'); xlabel('Fidelity');
mean_fidelity = mean(cell2mat(T_jitter_fidelity.fidelity));
std_fidelity = std(cell2mat(T_jitter_fidelity.fidelity));

text_str = sprintf("Mean \\pm STD = %0.2f \\pm %0.2f", mean_fidelity, std_fidelity);
text(0.4,0.4, text_str);
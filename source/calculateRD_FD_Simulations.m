% see the histograms.
close all;
close all;
clear all;
clear all;

d1 = load('histogramDataRes11.data');
d2 = load('histogramDataRes12.data');
d3 = load('histogramDataRes13.data');
d4 = load('histogramDataRes14.data');
d5 = load('histogramDataRes15.data');
d6 = load('histogramDataRes16.data');

% Normalisation with probability not pdf.
figure;
h1 = histogram(d1, 10000, 'Normalization','probability');
xlim([0 5]);
figure;
h2 = histogram(d2, 10000, 'Normalization','probability');
xlim([0 5]);
figure;
h3 = histogram(d3, 10000, 'Normalization','probability');
xlim([0 5]);
figure;
h4 = histogram(d4, 10000, 'Normalization','probability');
xlim([0 5]);
figure;
h5 = histogram(d5, 10000, 'Normalization','probability');
xlim([0 5]);
figure;
h6 = histogram(d6, 10000, 'Normalization','probability');
xlim([0 5]);

% make sure that you have area as 1. tjey are. now the covarnace which is RD.

d21 = load('histogramDataRes21.data');
d22 = load('histogramDataRes22.data');
d23 = load('histogramDataRes23.data');
d24 = load('histogramDataRes24.data');
d25 = load('histogramDataRes25.data');
d26 = load('histogramDataRes26.data');

% Normalisation with probability not pdf.
figure;
h21 = histogram(d21, 10000, 'Normalization','probability');
xlim([0 5]);
figure;
h22 = histogram(d22, 10000, 'Normalization','probability');
xlim([0 5]);
figure;
h23 = histogram(d23, 10000, 'Normalization','probability');
xlim([0 5]);
figure;
h24 = histogram(d24, 10000, 'Normalization','probability');
xlim([0 5]);
figure;
h25 = histogram(d25, 10000, 'Normalization','probability');
xlim([0 5]);
figure;
h26 = histogram(d26, 10000, 'Normalization','probability');
xlim([0 5]);


rd1 = std(h1.Values);
rd2 = std(h2.Values);
rd3 = std(h3.Values);
rd4 = std(h4.Values);
rd5 = std(h5.Values);
rd6 = std(h6.Values);

rd21 = std(h21.Values);
rd22 = std(h22.Values);
rd23 = std(h23.Values);
rd24 = std(h24.Values);
rd25 = std(h25.Values);
rd26 = std(h26.Values);

% fractals. rd1, rd2 are at m = 8 (1 mm resolution). rd21, rd22 are at m_ref = 1 (0.5 mm resolution)
fd1 = 1.0 - log(rd1/rd21) / log(8)
fd2 = 1.0 - log(rd2/rd22) / log(8)
fd3 = 1.0 - log(rd3/rd23) / log(8)
fd4 = 1.0 - log(rd4/rd24) / log(8)
fd5 = 1.0 - log(rd5/rd25) / log(8)
fd6 = 1.0 - log(rd6/rd26) / log(8)



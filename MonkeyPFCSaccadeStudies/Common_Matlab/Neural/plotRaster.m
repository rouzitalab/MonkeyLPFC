function plotRaster(cTrial)

[nSamples, nUnits] = size(cTrial.raster);
trTVec = (1:nSamples) - cTrial.eventTime;
raster = full(cTrial.raster);
raster = ((1:nUnits)' * ones(1, nSamples))' .* raster;
plot(trTVec, raster, 'b.', 'MarkerSize', 20),
hold on
plot([0 0], [0 nUnits], 'k--')
hold off
ylim([0.5 nUnits+0.1]);
xlim([min(trTVec) max(trTVec)])
axis ij
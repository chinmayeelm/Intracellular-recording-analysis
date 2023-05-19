rampData = readlines('rampDataList.txt');
dataDirectory = extractBefore(rampData, '_');
filename = extractAfter(rampData, "_");
nfiles = length(filename);
adaptation_coeff = zeros(1, nfiles);
rsquare = zeros(1, nfiles);

for irow = 1:nfiles
    
    P = getStructP(dataDirectory(irow), filename(irow));
    rampGCFRplots(P);
    [adaptation_coeff(irow), rsquare(irow)] = calcAdaptation(P);
    

end
function logTau = findSurfaceSpatial(CSF,FREQ)
% 
    logGain=CSF(:,1);
    logCenter=CSF(:,2);
    
    octaveWidth=10.^CSF(:,3);
    logTrunc=10.^CSF(:,4);
    logWidth = [octaveWidth.*log10(2)]./2;
    
    tauDecay = .5;
    K = log10(tauDecay);
    
    logP = logGain + K.*[(1./logWidth).*(FREQ - logCenter)].^2;

    truncLevel = logGain - logTrunc;

    leftCSF = [(logP < truncHalf) & (FREQ < logCenter)].*truncLevel;
    rightCSF = [(logP >= truncHalf) | (FREQ > logCenter)].*logP;

    logTau = -(leftCSF + rightCSF);

    logTau(find(logTau>0))=0;
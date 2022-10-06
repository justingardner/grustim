% find index for staircase
function idx = findCondIdx(staircaseTable,backLum,noiseLum,stimLum,stimDur,stimStd,stimColor)
    idx = (staircaseTable.backLum == backLum);
    idx = idx & (staircaseTable.noiseLum == noiseLum);
    idx = idx & (staircaseTable.stimLum == stimLum);
    idx = idx & (staircaseTable.stimDur == stimDur);
    idx = idx & (staircaseTable.stimStd == stimStd);
    idx = idx & (staircaseTable.stimColor == stimColor);
    idx = find(idx); % turn logical vector to index number
    
end
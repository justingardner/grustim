% find index for staircase
function idx = findCondIdx(staircaseTable,thistrial)
    idx = true(size(staircaseTable,1),1);
    T = staircaseTable;
    for pidx = 1:length(T.Properties.VariableNames)
        pname = T.Properties.VariableNames{pidx};
        if strcmp(pname, 'staircase')
            continue
        end
        if strcmp(pname, 'stimColor')
            idx  = idx & (T.(pname) == double(thistrial.(pname)));
        end
        
        if isfield(thistrial,pname)
            idx  = idx & (T.(pname) == thistrial.(pname));
        end
    end
    idx = find(idx); % turn logical vector to index number
end
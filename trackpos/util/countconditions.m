
function [changingfields, changingfield_vals] = countconditions(taskparams)
    changingfields = {};
    changingfield_vals = {};
    fieldnames = fields(taskparams);
    for fidx = 1:length(fieldnames)
        fieldname = fieldnames{fidx};
        if eval(['length(taskparams.' fieldname ') > 1'])
            changingfields{end+1} = fieldname;
            if strcmp(fieldname,'stimColor')
                changingfield_vals{end+1} = double(taskparams.stimColor);
            else
                changingfield_vals{end+1} = eval(['taskparams.' fieldname]);
            end
        end
    end
end
  
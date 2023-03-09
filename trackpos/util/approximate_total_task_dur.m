function totaldur = approximate_total_task_dur(task)
    % count trials
    numTrials        = 0; 
    totaltime        = 0;
    for phaseNum = 1:length(task{1})
        numTrials = numTrials + task{1}{phaseNum}.numTrials;
        totaltime = totaltime + task{1}{phaseNum}.numTrials * sum(task{1}{phaseNum}.segmax);
    end
    totaldur = totaltime/60/60;
    disp(['Approx task duration = ' num2str(totaldur) ' hours']);
end


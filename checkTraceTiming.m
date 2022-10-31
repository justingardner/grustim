%    checkTraceTiming.m
%
%       
%       usage: checkTraceTiming
%       by: josh wilson
%       date: October 2022
%       
%       Takes a stim file and plots the timing of the triggers
%      
%       Example Usage:
%           checkTraceTiming('stimFile=testTrigger')
%           

function checkTraceTiming(varargin)

% get stimfile name and load
getArgs(varargin);
load(stimFile)

%get the times
traces = makeTraces(myscreen) 
triggerTimes = myscreen.events.time(myscreen.events.tracenum == 1);
triggerTimes = triggerTimes - triggerTimes(1); %subtract timing of first trigger
triggerTimes = triggerTimes(2:end); %take out first backtick timing (experiment start)

% plot the trigger times
figure,
subplot(1,2,1); scatter(1:length(triggerTimes),triggerTimes);
hold on, plot([1:length(triggerTimes)], [triggerTimes(1):(triggerTimes(end)-triggerTimes(1))/(length(triggerTimes)-1):triggerTimes(end)]);
xlabel('Trigger number');ylabel('Time (s)');

% plot trigger times as difference from expect timing assuming all are equal
slope = (triggerTimes(end)-triggerTimes(1))/(length(triggerTimes)-1);

predictions = [triggerTimes(1):(triggerTimes(end)-triggerTimes(1))/(length(triggerTimes)-1):triggerTimes(end)];

subplot(1,2,2); scatter(1:length(triggerTimes), triggerTimes - predictions);

xlabel('Trigger number');ylabel('Difference from expected time (s)');
keyboard
% football_rodrigo_1.m
%
%      usage: football_rodrigo_1
%         by: justin gardner
% updated by: Rodrigo Sigala
%       date: 08/2010
%    purpose: program of distractor (saccade or flicker) contingent change blindness
%     input :
%             pathType: root path where all files are located
%               0 -'/Users/rodrigo/matlab/sigala/' (default)
%               1 -'~/proj/sigala/'
%
function football(varargin)

magnetRun = [];
mainPath = [];
getArgs(varargin,{'mainPath=~/grustim/','magnetRun=0'});


%%%%%%%%%%%%%%%%%%%%% Load necessary files %%%%%%%%%%%%%%%%%%%%%
%Image of a football (distractors)
footballImageFilename = fullfile(mainPath,'images/changeBlindness/ball.tif');

%Grid image (possibly used as distractor)
distractImageFilename = fullfile(mainPath,'images/changeBlindness/grid.tif');

%Invertee image of a football (necessaryo to create the flicker)
footballInvImageFilename = fullfile(mainPath,'images/changeBlindness/ball_inv.tif');

%File containing the positions of the distractors that increase the
%blindness effect
bckgConfigurationFilename = fullfile(mainPath,'football_config.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% initalize the screen
myscreen = initScreen;

task{1}.waitForBacktick = 1;
% There are 5 segments in the task
%{n}: task{1}.perceptualBlocker
% Segment 1, {0}  one football turns red and the subject has to saccade to that target
%            {1}  All footballs are displayed

% Segment 2, {0} same as segment 1
%            {1} Some balls start to flicker and balls start to flicker


% Segment 3, {0} the second football turns red to cue the location of the target
%            {1} Same as Segment 2

% Segment 4, {0} the added target turns on if this is not a distractorContingent trial
%            {0} the added target continues on without flickering
%            distractor
% Segment 5, is when a single football comes out in green or red to tell the subject
% whether they got the trial correct or not.
task{1}.segmin = [1 1 1 1 1];
if magnetRun
  task{1}.segmax = [1 1 1 1 11];
  task{1}.synchToVol = [0 0 0 0 1];
else
  task{1}.segmax = [1 1 1 1 1];
  task{1}.synchToVol = [0 0 0 0 0];
end  
%task{1}.segquant = [0 0.02 0 0 0];

% we allow the observer to respond only in segment 3-5
task{1}.getResponse = [0 1 1 1 1];

%%%%%%%%%%% Randomized parameters:

% distractorContingent controls whether the change will occur during a
% distractor  or not
task{1}.randVars.uniform.distractorContingent = [0 1 1 1];

% presentAbsent controls whether the target will appear or not.
task{1}.randVars.uniform.presentAbsent = [0 1];

%perceptualBlocker controls the source of the change blindness 0: Saccade;
%1: Flickering balls;
task{1}.randVars.uniform.perceptualBlocker = [1];


%The following parameters are only relevant if perceptualBlocker is 1
%(flickering)

%useFixedBckgs: Determines whether the distractors will be defined randomly
%(0) or based in a set of predefined templates previously stored (>0)
%1: Templates will be used for the positions of ALL distractors
%2: Templates will be used to define positions of the balls in the same
%quadrant of the target (to keep "target quadrants" the same)
task{1}.randVars.uniform.useFixedBckgs = [2];


%The following parameters are only relevant if useFixedBackds
%is set to 0 (not templates used)

% distractorPlacement controls whether the flicker stimulis is restricked
%to the contralateral side of the hemifield
%distractorPlacement: 0: Totally random
%                     1: Placed on the half of the visual field contrary to
%                     the position of the target.
%                     2: Placed on the other .3 quadrants relative to the
%                     position of the target.
task{1}.randVars.uniform.distractorPlacement = [0];

%distributeDistractors: try to locate the distractors (And flickering balls
%in different locations of the visual field to increase the  effect of
%cluter
task{1}.randVars.uniform.distributeDistractors = [1];

%sequenceFlickers determines whether the flickerng balls appear in sequence
%or all together
task{1}.randVars.uniform.sequenceFlickers = [0];

%Gives feedback to the subject about their performance
task{1}.randVars.uniform.giveFeedback = [0];

% random makes the parameters above show in a block randomized order.
task{1}.random = 1;


% we run 100 trials and then stop
task{1}.numTrials = 100;


% initialize the task
for phaseNum = 1:length(task)
    [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@drawStimulusCallback,@responseCallback);
end

% init the stimulus.
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

% Call myInitStimulus which just loads the image of the footbal 221l and
% creates textures with different color footballs.
stimulus = myInitStimulus(stimulus,myscreen,footballImageFilename,footballInvImageFilename,distractImageFilename, bckgConfigurationFilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;

flickerReleased=1;
flickerCounter=1;

while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    % task{1}
    %  task{1}.thistrial
    % update the task
    [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
    
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

% quick analysis of the data. Only run this if we have enough trials.
if task{1}.trialnum > 10
    % get the task parameters
    e = getTaskParameters(myscreen,task);
    % get the value of distractorContingent and presentAbsent for each trial
    distractorContingent = e.randVars.distractorContingent;
    presentAbsent = e.randVars.presentAbsent;
    % get the response on each trial
    response = (e.response == 1);
    % compute distractorContingent hits and falseAlarams
    hits = sum(presentAbsent(find(distractorContingent)) & response(find(distractorContingent)));
    hitsn = sum(presentAbsent(find(distractorContingent)));
    falseAlarms = sum(~presentAbsent(find(distractorContingent)) & response(find(distractorContingent)));
    falseAlarmsN = sum(~presentAbsent(find(distractorContingent)));
    % display hits and false alarms
    disp(sprintf('(football) distractorContingent hits: %0.2f%% (%i/%i) false alarms: %0.2f%% (%i/%i)',100*hits/hitsn,hits,hitsn,100*falseAlarms/falseAlarmsN,falseAlarms,falseAlarmsN));
    % same for non distractorContingent
    hits = sum(presentAbsent(find(~distractorContingent)) & response(find(~distractorContingent)));
    hitsn = sum(presentAbsent(find(~distractorContingent)));
    falseAlarms = sum(~presentAbsent(find(~distractorContingent)) & response(find(~distractorContingent)));
    falseAlarmsN = sum(~presentAbsent(find(~distractorContingent)));
    disp(sprintf('(football) Not distractorContingent hits: %0.2f%% (%i/%i) false alarms: %0.2f%% (%i/%i)',100*hits/hitsn,hits,hitsn,100*falseAlarms/falseAlarmsN,falseAlarms,falseAlarmsN));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

if task.thistrial.thisseg == 1
    
    % after the fifth trial, we stop showing where the eye position is and
    % changing the color of the football when the subject has fixated the target.
    if task.trialnum > 5
        stimulus.displayEyeContingent = 0;
    end
    
    % variable just controls whether a target should be added. In the first segment
    % we always omit the target - unless, presentAbsent is set to -1 in which
    % case we will disappear rather than appear the target
    if task.thistrial.presentAbsent == -1
      stimulus.addTarget = 1;
    else
      stimulus.addTarget = 0;
    end
    
    % clear screen to gray
    mglClearScreen(0.5);
    
    %targetOptionsX  = [-12 -6 6 12];
    %targetOptionsY  = [-8 -4 4 8];
    
    targetOptionsX  = [-12  12];
    targetOptionsY  = [-8   8];
    
    targetSize = 4.5;
    
    %Controls the global display properties
    %[a b c] : a: shift in the horizontal axis (Degs); b: shift in the vertical
    %axis (Degs); c: scale factor (0:1)
    stimulus.displayParam = [0 0 1];
    
    % get a set of random stimulus locations and sizes
    % Also construct a vector of textures (distractors)
    stimulus.textVect=[];
    stimulus.invTextVect=[];
    
    %Use predined backgrounds
    if task.thistrial.useFixedBckgs==1
        disp(['Backgnd set: ' num2str(stimulus.bckgndConf)]);
        stimulus.loc = squeeze(stimulus.confBackgs(stimulus.bckgndConf,:,:));
        
        for i=1:stimulus.nDistract
            stimulus.textVect = [stimulus.textVect stimulus.tex];
            stimulus.invTextVect = [stimulus.invTextVect stimulus.invtex];
        end
    elseif task.thistrial.useFixedBckgs==2
        
        stimulus.bckgndConf = 9;  %Configuration defined arbitrarily because it works good
        ballsInQuad=intersect(find (stimulus.confBackgs(stimulus.bckgndConf,:,1)<0),find (stimulus.confBackgs(stimulus.bckgndConf,:,2)>0));
        stimulus.loc = squeeze(stimulus.confBackgs(stimulus.bckgndConf,:,:));
        
        for i = 1:stimulus.nStimuli
            
            validated = 0;  %Avoids a seen that has no overlapp with the target stimulus
            sameQuad=0;
            
            if i==1
                if length(targetOptionsX)==4
                    stimulus.loc(i,1:4) = [targetOptionsX(round(mod(rand/.03,3))+1) targetOptionsY(round(mod(rand/.03,3))+1) targetSize targetSize];
                else
                    stimulus.loc(i,1:4) = [targetOptionsX(round(rand)+1) targetOptionsY(round(rand)+1) targetSize targetSize];
                end
            else
                while ~validated
                    
                    s = rand*4+2;
                    
                    if (stimulus.loc(1,1)<0 && stimulus.loc(1,2)>0)  %First quad
                        if ~isempty(intersect(i,ballsInQuad))  %if is already defined
                            %correctVector=[-1*sign(squeeze(stimulus.confBackgs(stimulus.bckgndConf,i,1:2))); 1; 1];
                            stimulus.loc(i,1:4) = squeeze(stimulus.confBackgs(stimulus.bckgndConf,i,:));
                            validated=1;
                        else
                            stimulus.loc(i,1:4) = [rand*stimulus.deviceWidth-stimulus.deviceWidth/2 rand*stimulus.deviceHeight-stimulus.deviceHeight/2 s s];
                            if xor(stimulus.loc(1,1)>0,stimulus.loc(i,1)>0) || xor(stimulus.loc(1,2)>0,stimulus.loc(i,2)>0)
                                validated=1;
                            else
                                validated=0;
                            end
                        end
                        
                    elseif (stimulus.loc(1,1)>0 && stimulus.loc(1,2)>0)  %Sec quad
                        if ~isempty(intersect(i,ballsInQuad))  %if is already defined
                            correctVector=[sign(squeeze(stimulus.confBackgs(stimulus.bckgndConf,i,1:2))); 1; 1];
                            stimulus.loc(i,1:4) = squeeze(stimulus.confBackgs(stimulus.bckgndConf,i,:)).*correctVector;
                            validated=1;
                        else
                            stimulus.loc(i,1:4) = [rand*stimulus.deviceWidth-stimulus.deviceWidth/2 rand*stimulus.deviceHeight-stimulus.deviceHeight/2 s s];
                            if xor(stimulus.loc(1,1)>0,stimulus.loc(i,1)>0) || xor(stimulus.loc(1,2)>0,stimulus.loc(i,2)>0)
                                validated=1;
                            else
                                validated=0;
                            end
                        end
                        
                    elseif (stimulus.loc(1,1)>0 && stimulus.loc(1,2)<0)  %Third quad
                        if ~isempty(intersect(i,ballsInQuad))  %if is already defined
                            correctVector=[-1*sign(squeeze(stimulus.confBackgs(stimulus.bckgndConf,i,1:2))); 1; 1];
                            stimulus.loc(i,1:4) = abs(squeeze(stimulus.confBackgs(stimulus.bckgndConf,i,:))).*correctVector;
                            validated=1;
                        else
                            stimulus.loc(i,1:4) = [rand*stimulus.deviceWidth-stimulus.deviceWidth/2 rand*stimulus.deviceHeight-stimulus.deviceHeight/2 s s];
                            if xor(stimulus.loc(1,1)>0,stimulus.loc(i,1)>0) || xor(stimulus.loc(1,2)>0,stimulus.loc(i,2)>0)
                                validated=1;
                            else
                                validated=0;
                            end
                        end
                    else %Fourth quad
                        if ~isempty(intersect(i,ballsInQuad))  %if is already defined
                            correctVector=[-1*sign(squeeze(stimulus.confBackgs(stimulus.bckgndConf,i,1:2))); 1; 1];
                            stimulus.loc(i,1:4) = squeeze(stimulus.confBackgs(stimulus.bckgndConf,i,:)).*correctVector;
                            validated=1;
                        else
                            stimulus.loc(i,1:4) = [rand*stimulus.deviceWidth-stimulus.deviceWidth/2 rand*stimulus.deviceHeight-stimulus.deviceHeight/2 s s];
                            if xor(stimulus.loc(1,1)>0,stimulus.loc(i,1)>0) || xor(stimulus.loc(1,2)>0,stimulus.loc(i,2)>0)
                                validated=1;
                            else
                                validated=0;
                            end
                        end
                    end
                    
                end
            end
            
            
            %Is a flickering ball
            if i>stimulus.nStimuli-stimulus.nDistract
                stimulus.textVect = [stimulus.textVect stimulus.tex];
                stimulus.invTextVect = [stimulus.invTextVect stimulus.invtex];
            end
        end
        
        stimulus.loc(:,1:2) = stimulus.loc(:,1:2) * stimulus.displayParam(3);  %Scale screen if necessary
        stimulus.loc(:,1) = stimulus.loc(:,1) + stimulus.displayParam(1);
        stimulus.loc(:,2) = stimulus.loc(:,2) + stimulus.displayParam(2);
        
        
    else
        for i = 1:stimulus.nStimuli
            
            validated = 0;
            sameQuad=0;
            
            if i==1
                if length(targetOptionsX)==4
                    stimulus.loc(i,1:4) = [targetOptionsX(round(mod(rand/.03,3))+1) targetOptionsY(round(mod(rand/.03,3))+1) targetSize targetSize];
                else
                    stimulus.loc(i,1:4) = [targetOptionsX(round(rand)+1) targetOptionsY(round(rand)+1) targetSize targetSize];
                end
            else
                while ~validated
                    s = rand*4+2;
                    
                    if task.thistrial.distractorPlacement==1 && i>stimulus.nStimuli-stimulus.nDistract
                        if stimulus.loc(1,1)>0  %right hemifield
                            %generate flickering distractors on the LEFT side
                            stimulus.loc(i,1:4) = [-abs(rand*stimulus.deviceWidth-stimulus.deviceWidth/2) rand*stimulus.deviceHeight-stimulus.deviceHeight/2 s s];
                            
                        else %left hemifield
                            %generate flickering distractors on the RIGHT side
                            stimulus.loc(i,1:4) = [abs(rand*stimulus.deviceWidth-stimulus.deviceWidth/2) rand*stimulus.deviceHeight-stimulus.deviceHeight/2 s s];
                        end
                        
                    elseif task.thistrial.distractorPlacement==2 && i>stimulus.nStimuli-stimulus.nDistract
                        stimulus.loc(i,1:4) = [rand*stimulus.deviceWidth-stimulus.deviceWidth/2 rand*stimulus.deviceHeight-stimulus.deviceHeight/2 s s];
                        %Check if distractor is in the same quadrant but is not
                        %allowd
                        if xor(stimulus.loc(1,1)>0,stimulus.loc(i,1)>0) || xor(stimulus.loc(1,2)>0,stimulus.loc(i,2)>0)
                            sameQuad=0;
                        else
                            sameQuad=1;
                        end
                    elseif task.thistrial.distributeDistractors
                        switch mod(i,4)+1
                            case 1,
                                signX=-1;
                                signY=1;
                            case 2,
                                signX=1;
                                signY=1;
                            case 3,
                                signX=1;
                                signY=-1;
                            case 4,
                                signX=-1;
                                signY=-1;
                        end
                        stimulus.loc(i,1:4) = [signX*abs(rand*stimulus.deviceWidth-stimulus.deviceWidth/2) signY*abs(rand*stimulus.deviceHeight-stimulus.deviceHeight/2) s s];
                    else
                        stimulus.loc(i,1:4) = [rand*stimulus.deviceWidth-stimulus.deviceWidth/2 rand*stimulus.deviceHeight-stimulus.deviceHeight/2 s s];
                    end
                    
                    %                 [stimulus.loc(1,1:2);stimulus.loc(i,1:2)]
                    %                 (abs(stimulus.loc(i,1)-stimulus.loc(1,1)) + abs(stimulus.loc(i,2)-stimulus.loc(1,2)))
                    
                    %If target and distractor overlap
                    if (abs(stimulus.loc(i,1)-stimulus.loc(1,1)) + abs(stimulus.loc(i,2)-stimulus.loc(1,2))) < 8
                        validated=0;
                    else
                        validated=1;
                        if sameQuad
                            validated=0;
                        end
                    end
                end
                
                %Is a flickering ball
                if i>stimulus.nStimuli-stimulus.nDistract
                    stimulus.textVect = [stimulus.textVect stimulus.tex];
                    stimulus.invTextVect = [stimulus.invTextVect stimulus.invtex];
                end
            end
        end
    end
    
    % now get the location of two of the targets. The subject will be cued to saccade
    % to the first target loc and then the second tar qget loc - i.e. the football
    % at targetLoc1 will change red and then the football at targetLoc2 will change to red.
    stimulus.targetLoc1 = stimulus.loc(stimulus.nStimuli,:);
    stimulus.targetLoc2 = stimulus.loc(stimulus.nStimuli-1,:);
    % the current targetLoc for the first segment will be targetLoc1
    stimulus.targetLoc = stimulus.targetLoc1;
    
    % This is a bit redundant, but keep a cell array of all the target locations
    % for each trial - this is so that later we could do an analysis based
    % on the exact location and sizes of everythin gthat was shown.
    stimulus.locs{task.trialnum}.targetLoc1 = stimulus.targetLoc1;
    stimulus.locs{task.trialnum}.targetLoc2 = stimulus.targetLoc2;
    stimulus.locs{task.trialnum}.newTarget = stimulus.loc(1,:);
    stimulus.locs{task.trialnum}.locs = stimulus.loc;
elseif task.thistrial.thisseg == 3
    % on segment 3 we change the location of the saccade target (turn red) to
    % the second target, i.e. targetLoc2
    stimulus.targetLoc = stimulus.targetLoc2;
end

if task.thistrial.perceptualBlocker
  % here, we flip the sign of stimulus.addTarget. That way
  % if we want to disappear the target, we just need to set
  % addTarget to 1 in the beginning. If we want to appear
  % the target we set addTarget to 0 in the beginning
  if task.thistrial.presentAbsent 
    if  task.thistrial.distractorContingent
      if task.thistrial.thisseg == 2 && stimulus.thistrial.flickerNumber>1
	stimulus.addTarget = xor(stimulus.addTarget,1);
	stimulus.thistrial.tFlick=mglGetSecs(stimulus.thistrial.t0Flick);
      end
    elseif task.thistrial.thisseg == 3
      stimulus.addTarget = xor(stimulus.addTarget,1);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

% we check the response for correct or incorrect
if task.thistrial.gotResponse < 1
    if task.thistrial.presentAbsent && (task.thistrial.whichButton == 1)
        % this is a hit
        stimulus.type(task.trialnum) = 'H';
        stimulus.correct(task.trialnum) = 1;
        
    elseif task.thistrial.presentAbsent && (task.thistrial.whichButton == 2)
        % this is a miss
        stimulus.type(task.trialnum) = 'M';
        stimulus.correct(task.trialnum) = 0;
    elseif ~task.thistrial.presentAbsent && (task.thistrial.whichButton == 2)
        % this is a correct reject
        stimulus.type(task.trialnum) = 'C';
        stimulus.correct(task.trialnum) = 1;
    elseif ~task.thistrial.presentAbsent && (task.thistrial.whichButton == 1)
        % this is a false alarm
        stimulus.type(task.trialnum) = 'F';
        stimulus.correct(task.trialnum) = 0;
    end
end

% display to text windown what happened
disp(sprintf('Trial %i: (%s) distractorContingent: %i',task.trialnum,stimulus.type(task.trialnum),task.thistrial.distractorContingent));

e = getTaskParameters(myscreen,task);


%Add trial type
response = e.response(e.nTrials)==1;
distractorContingent = e.randVars.distractorContingent(e.nTrials);
presentAbsent = e.randVars.presentAbsent(e.nTrials);

% Classify cuadrant of the target stimulus:
% 2_|_1
% 4_|_3
xVal = (stimulus.loc(1,1)<0)+1;
yVal = (stimulus.loc(1,2)<0)+1;
cuadrant = mod(3^xVal + 2^yVal,5)+1;
trialCond = 4; %Conditions == number of quadrants

%[xVal yVal cuadrant]

if distractorContingent
    if presentAbsent&&response %Hit`2
        trialType = cuadrant;
    elseif presentAbsent&&~response  %Mistaken
        trialType = trialCond + cuadrant;
    elseif ~presentAbsent&&response %False alarm
        trialType = trialCond*2 + cuadrant;
    else  %Correct
        trialType = trialCond*3 + cuadrant;
    end
else
    if presentAbsent&&response %Hits
        trialType = trialCond*4 + cuadrant;
    elseif presentAbsent&&~response  %Mistaken
        trialType = trialCond*5 + cuadrant;
    elseif ~presentAbsent&&response %False alarm
        trialType = trialCond*6 + cuadrant;
    else  %Correct
        trialType = trialCond*7 + cuadrant;
    end
end

%task.thistrial.quadrantTarget = cuadrant;
%task.thistrial.trialTypeAll = trialType;
task.randVars.quadrantTarget(e.nTrials) = cuadrant;
task.randVars.trialTypeAll(e.nTrials) = trialType;



%disp(['Timing of target on after the beggining of flickering: ' num2str(e.trials(e.nTrials).segtime(3)-e.trials(e.nTrials).segtime(2))])
%disp(['distractorPlacement:' num2str(task.thistrial.distractorPlacement) ])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = drawStimulusCallback(task, myscreen)

global stimulus

% this variable just controls how close in degrees to the edge of the football
% the saccade has to be to considered on target.
fixAccuracy = 5;

% clear the screen
mglClearScreen;

%fixationCross size
fixationCrossSize = 1;


% on the fifth segment, we are just going to show a green or red football
% depending on whether the subject got the trial correct or not.
if task.thistrial.thisseg == 5
    
    if task.thistrial.gotResponse
        mglFixationCross(fixationCrossSize,1.5,[1 0 0],stimulus.displayParam(1:2));
    else
        %mglFixationCross;
        mglFixationCross(fixationCrossSize,1.5,[1 1 1],stimulus.displayParam(1:2));
    end
    %return
    if task.thistrial.giveFeedback
        if (length(stimulus.correct) >= task.trialnum)
            if stimulus.correct(task.trialnum)
                mglBltTexture(stimulus.greentex,[0 0]) ;
            else
                mglBltTexture(stimulus.redtex,[0 0]);
            end
        end
    else
        return
    end
    %return
end

% on every other semgent show all of the footballs except the first one
% (target stimulus)
for i = 2:stimulus.nStimuli
    % blt the football texture to the correct location and size
    mglBltTexture(stimulus.tex,stimulus.loc(i,:));
end

if task.thistrial.perceptualBlocker==0
    
    % check to see if the subject has fixated targetLoc1
    if ((myscreen.eyetracker.eyepos(1) > (stimulus.targetLoc1(1)-stimulus.targetLoc1(3)/2-fixAccuracy/2)) &&...
            (myscreen.eyetracker.eyepos(1) < (stimulus.targetLoc1(1)+stimulus.targetLoc1(3)/2+fixAccuracy/2)) &&...
            (myscreen.eyetracker.eyepos(2) > (stimulus.targetLoc1(2)-stimulus.targetLoc1(4)/2-fixAccuracy/2)) &&...
            (myscreen.eyetracker.eyepos(2) < (stimulus.targetLoc1(2)+stimulus.targetLoc1(4)/2+fixAccuracy/2)))
        stimulus.targetFixed1 = 1;
    else
        stimulus.targetFixed1 = 0;
    end
    
    % check to see if the subject has fixated targetLoc
    if ((myscreen.eyetracker.eyepos(1) > (stimulus.targetLoc2(1)-stimulus.targetLoc2(3)/2-fixAccuracy/2)) &&...
            (myscreen.eyetracker.eyepos(1) < (stimulus.targetLoc2(1)+stimulus.targetLoc2(3)/2+fixAccuracy/2)) &&...
            (myscreen.eyetracker.eyepos(2) > (stimulus.targetLoc2(2)-stimulus.targetLoc2(4)/2-fixAccuracy/2)) &&...
            (myscreen.eyetracker.eyepos(2) < (stimulus.targetLoc2(2)+stimulus.targetLoc2(4)/2+fixAccuracy/2)))
        stimulus.targetFixed2 = 1;
    else
        stimulus.targetFixed2 = 0;
    end
    
    % display the saccade target in red.
    mglBltTexture(stimulus.redtex,stimulus.targetLoc);
    
    % this code displays the fixation accuracy. It turns the footballs green
    % when the subject has fixaated them and shows a little black spot
    % for where the eye position is. After five trials this initSegmentCallback
    % turns off displayEyeContingent so that this doesn't display anymore
    
    if stimulus.displayEyeContingent
        % display in green if we have saccaded to targetLoc1
        if stimulus.targetFixed1
            mglBltTexture(stimulus.greentex,stimulus.targetLoc1);
        end
        
        % display in green i f we have saccaded to targetLoc2
        if (task.thistrial.thisseg >= 2) && stimulus.targetFixed2
            %if (task.thistrial.thisseg >= 2)
            mglBltTexture(stimulus.greentex,stimulus.targetLoc2);
            %mglBltTexture(stimulus.invtex,stimulus.loc(stimulus.nStimuli,:))
        end
        
        % display eye position
        mglGluDisk(myscreen.eyetracker.eyepos(1),myscreen.eyetracker.eyepos(2), 0.1,[0.3 05 0.8]);
    end
    % if this is an eyeContingent trial,
    % now, if we are supposed to add a target
    if task.thistrial.presentAbsent
        % if this is a distractor contingent trial then we display the added target
        % when the subject is no longer fixating target1 (i.e. it turns on
        % as the subject moves their eyes away from targe tLoc1).
        if task.thistrial.distractorContingent
            if (task.thistrial.thisseg == 3) && (stimulus.targetFixed1 == 0)
                stimulus.addTarget = 1;
            end
            % otherwise just display the target in segment4.
        elseif task.thistrial.thisseg == 4
            stimulus.addTarget = 1;
        end
    end
    
    
    %%%%%%%%%%%%% Flickering stimulus
elseif task.thistrial.perceptualBlocker>0
    
    %stimulus.thistrial.flickerDuration = .1;
    stimulus.thistrial.flickerDuration = 6;  %Defined in screen refresh cycles
    
    if (task.thistrial.thisseg == 1)
        stimulus.thistrial.t0Flick=mglGetSecs;
        stimulus.thistrial.flickerNumber = 0;
    end
    
    
    %Adds a distractor
    if (task.thistrial.thisseg == 2) && task.thistrial.distractorContingent &&  stimulus.thistrial.flickerNumber<stimulus.thistrial.flickerDuration
        if task.thistrial.perceptualBlocker<2  %Flickering balls restricted in position
            if task.thistrial.sequenceFlickers %Sequence of balls (Section not revised)
                flickBall=1;
                stimulus.thistrial.tFlick=mglGetSecs(stimulus.thistrial.t0Flick);
                while flickBall<stimulus.nDistract && stimulus.thistrial.tFlick<stimulus.thistrial.flickerDuration  %Flickering
                    mglBltTexture(stimulus.invtex,stimulus.loc(stimulus.nStimuli-stimulus.nDistract+flickBall,:));
                    mglFlush;
                    mglBltTexture(stimulus.tex,stimulus.loc(stimulus.nStimuli-stimulus.nDistract+flickBall,:));
                    mglFlush;
                    flickBall=flickBall+1;
                end
            else
                if mod(stimulus.thistrial.flickerNumber,2)
                    mglBltTexture(stimulus.invTextVect,stimulus.loc(stimulus.nStimuli-stimulus.nDistract+1:stimulus.nStimuli,:));
                else
                    mglBltTexture(stimulus.textVect,stimulus.loc(stimulus.nStimuli-stimulus.nDistract+1:stimulus.nStimuli,:));
                end
                stimulus.thistrial.flickerNumber=stimulus.thistrial.flickerNumber+1;
            end
        else  % Flickering quads
        end
    end
    
    %if ((task.thistrial.thisseg == 3)) && task.thistrial.distractorContingent && stimulus.thistrial.tFlick<stimulus.thistrial.flickerDuration  %Flickering
    if ((task.thistrial.thisseg == 3)) && task.thistrial.distractorContingent && stimulus.thistrial.flickerNumber<stimulus.thistrial.flickerDuration  %Flickering
        if mod(stimulus.thistrial.flickerNumber,2)
            mglBltTexture(stimulus.invTextVect,stimulus.loc(stimulus.nStimuli-stimulus.nDistract+1:stimulus.nStimuli,:));
        else
            mglBltTexture(stimulus.textVect,stimulus.loc(stimulus.nStimuli-stimulus.nDistract+1:stimulus.nStimuli,:));
        end
        stimulus.thistrial.flickerNumber=stimulus.thistrial.flickerNumber+1;
    end
    
end


% draw the added target
if stimulus.addTarget
    mglBltTexture(stimulus.tex,stimulus.loc(1,:));
end


if task.thistrial.gotResponse
    mglFixationCross(fixationCrossSize,1.5,[1 0 0 ],stimulus.displayParam(1:2));
else
    %mglFixationCross;
    mglFixationCross(fixationCrossSize,1.5,[1 1 1],stimulus.displayParam(1:2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,footballImageFilename,footballInvImageFilename,distractImageFilename,bckgConfigurationFilename)

% load the stimulus
im1 = imread(footballImageFilename);
im2(:,:,1) = double(flipud(im1(:,:,1)));
im2(:,:,2) = im2(:,:,1);
im2(:,:,3) = im2(:,:,1);
im2(:,:,4) = im1(:,:,2);
stimulus.tex = mglCreateTexture(im2);

% make red one
im2(:,:,1) = double(flipud(im1(:,:,1)));
im2(:,:,2) = round(0.2*im2(:,:,1));
im2(:,:,3) = round(0.2*im2(:,:,1));
im2(:,:,4) = im1(:,:,2);
stimulus.redtex = mglCreateTexture(im2);

%myscreen

% make green one
im2(:,:,2) = double(flipud(im1(:,:,1)));
im2(:,:,1) = round(0.2*im2(:,:,2));
im2(:,:,3) = round(0.2*im2(:,:,2));
im2(:,:,4) = im1(:,:,2);
stimulus.greentex = mglCreateTexture(im2);

% grid
% im3 = imread(distractImageFilename);
% im4(:,:,1) = double(flipud(im3(:,:,1)));
% im4(:,:,2) = im4(:,:,1);
% im4(:,:,3) = im4(:,:,1);
% im4(:,:,4) = im3 (:,:,2);
% stimulus.disttex = mglCreateTextur    e(im4);

% inverted ball
im1_inv = imread(footballInvImageFilename);
im5(:,:,1) = double(flipud(im1_inv(:,:,1)));
im5(:,:,2) = im5(:,:,1);
im5(:,:,3) = im5(:,:,1);
im5(:,:,4) = im5 (:,:,2);
stimulus.invtex = mglCreateTexture(im5);


% size of screen (make a bit smaller than actual size to avoid
% having the sitmuli show off the edge of screen.
stimulus.deviceWidth = mglGetParam('deviceWidth')*3/4;
stimulus.deviceHeight = mglGetParam('deviceHeight')*3/4;

% number of stimuli
stimulus.nStimuli = 30;
stimulus.nDistract = round(stimulus.nStimuli*.4);

% will display the eye position in the beginning few trials
% but after 5 trials trun it off.
stimulus.displayEyeContingent = 1;

% default
stimulus.correct = [];

%Load some backgound configurations that work in confusing the subjects
%load('confs.mat','confBackgs');
%load('confs2.mat','confBackgs');
%load('~/proj/sigala/confs3.mat','confBackgs');
if isfile(bckgConfigurationFilename)
  disp(sprintf('(football_rodrigo_1) Loading configuration file %s',bckgConfigurationFilename));
  load(bckgConfigurationFilename,'confBackgs');
else
  disp(sprintf('(football_rodrigo_1) **** Could not find configuration file %s ****',bckgConfigurationFilename));
  keyboard
end

iniRow=16;
finRow=16;

%Modifying the backgnd
confBackgs(9,iniRow:finRow,:)=repmat([-9 4 5 5],finRow-iniRow+1,1);

stimulus.confBackgs = confBackgs;
%stimulus.goodConfs=[1 2 4 9 10]  %for confs1
%stimulus.goodConfs=[1 1 1 1 1 1]  %for confs1


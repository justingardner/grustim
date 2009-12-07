% readWords 
%
%      usage: myscreen = readWords(myscreen, filename)
%         by: justin gardner
%       date: 06/29/2006
%    purpose: read word files for fix task
%
%       e.g.:
%      
function myscreen = readWords(myscreen,filename)

global MGL

% open the file
fwords = fopen(filename,'r');
if (fwords == 0)
  disp(sprintf('UHOH: Could not open %s',filename));
  return
end
disp(sprintf('Using wordlist %s',filename));

% read all the words
[words nwords] = fscanf(fwords,'%s %i');

% close the file
fclose(fwords);

% if this is not divisible by six then we don't
% have a list of 5 letter words with single
% numbers, so complain
if (6*nwords/2 ~= size(words,1))
  disp(sprintf('UHOH: %s is not a list of 5 letter words',filename));
  return
end

%reshape the matrix appropriately
words = reshape(words,6,nwords/2)';

% get the words
myscreen.text.words = char(words);
myscreen.text.syllables = words(:,6);
myscreen.text.n = size(words,1);

disppercent(-inf,'Creating word textures');

% the mgl way
% * set defaults for fonts
% * create textures
mglTextSet(myscreen.text.font, myscreen.text.size,[1 1 1])

for wordnum = 1:myscreen.text.n
  disppercent(wordnum/myscreen.text.n);
  for casenum = 1:2
    % get the word we will display
    if (casenum == 2)
      thisword = upper(myscreen.text.words(wordnum,:));
    else
      thisword = lower(myscreen.text.words(wordnum,:));
    end
    
    % literally only need to do this
    myscreen.text.tex{wordnum,casenum} =  mglText(thisword);
  end
end
disppercent(inf);

% set up mode
myscreen.text.mode = -1;

% generate a couple random dot seeds
myscreen.text.randseed(1) = sum(100*clock);
myscreen.text.randseed(2) = sum(50*clock);

% set up task
myscreen.text.task = -1;
myscreen.text.this.gotresponse = 0;
myscreen = initTextTask(myscreen,1);

return

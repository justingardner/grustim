% mglFixationCross.m
%
%        $Id$
%      usage: mglFixationCross([width], [linewidth], [color], [origin]);
%             mglFixationCross( params ) 
%               params(1)=width
%               params(2)=linewidth
%               params(3:5)=color
%               params(6:7)=origin
%         by: jonas larsson
%       date: 
%  copyright: (c) 2006 Justin Gardner, Jonas Larsson (GPL see mgl/COPYING)
%    purpose: draws a fixation cross
%             with no arguments, draws a fixation cross at origin
%             (default width 0.2 with linewidth 1 in white at [0,0])
%       e.g.:
%
%mglOpen;
%mglVisualAngleCoordinates(57,[16 12]);
%mglFixationCross;
%mglFlush;
function mglFixationCrossArms( width, linewidth, color, origin );

if (exist('width','var') & length(width)>1)
  params=width;
  width=params(1);
  linewidth=params(2);
  if (length(params)>2)
    color=params(3);
    if (length(params)>3)
      origin=params(4:5);
    end
  end
end

if (~exist('width','var'))
  width=0.5;
end
if (~exist('origin','var'))
  origin=[0,0];
end
if (~exist('linewidth','var'))
  linewidth=1;
end
if (~exist('color','var'))
  color{1}=[0 0 0];
  color{2}=[0 0 0];
  color{3}=[0 0 0];
  color{4}=[0 0 0];
end

% up line
x0=origin(1);
x1=origin(1);
y0=origin(2);
y1=origin(2)+width;
mglLines2( x0,y0,x1,y1, linewidth, color{1} );
% down line
x0=origin(1);
x1=origin(1);
y0=origin(2);
y1=origin(2)-width;
mglLines2( x0,y0,x1,y1, linewidth, color{2} );
% left line
x0=origin(1);
x1=origin(1)-width;
y0=origin(2);
y1=origin(2);
mglLines2( x0,y0,x1,y1, linewidth, color{3} );
% right line
x0=origin(1);
x1=origin(1)+width;
y0=origin(2);
y1=origin(2);
mglLines2( x0,y0,x1,y1, linewidth, color{4} );



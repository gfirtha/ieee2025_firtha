function progbar(start,stop,value,varargin)
%PROGBAR   Show progress bar
%
%   PROGBAR(start,stop,value,'Key','Value',...) shows a progress bar in the
%   console window.
%
%   start    Value corresponding to an empty progress bar.
%   stop     Value corresponding to a full progress bar.
%   value    Current value of the progress bar.  On the first call of this
%            function, 'value' has to equal 'start' in order to trigger the
%            initialization of the progress bar, while on the last call it has
%            to equal 'stop' to trigger the erasing of the progress bar.
%   'Keep'   Keep displaying the full progress bar on completion.
%            'yes' or 'no'.
%            Default: 'no'.
%   'Width'  Total width (in characters) of the progress bar.
%            Default: console width with a maximum of 80 characters.

% Mattias Schevenels
% October 2004
% last modified: Peter Fiala 14.07.2014

% DECLARE PERSISTENT VARIABLES
persistent oldpos;
persistent tstart;  % start time

if stop <= start
    return;
end

% OBTAIN CONSOLE WIDTH
commandwindowsize=get(0,'CommandWindowSize');

% INTERPRETE KEY OPTIONS
keep=key('Keep','no',varargin{:});
width=key('Width',min(commandwindowsize(1),80),varargin{:});

% DETERMINE (NEW) POSITION
pos=round((value-start)/(stop-start)*(width-2));

% EXIT IF THE POSITION DID NOT CHANGE W.R.T. THE PREVIOUS CALL
% UNLESS THIS IS THE FIRST OR THE LAST CALL FOR THE CURRENT PROGRESS BAR
if (value~=stop) && (value~=start) && (pos==oldpos)
    return;
end

% ERASE OF THE PROGRESS BAR (IF REQUESTED) IF THIS IS THE LAST CALL
% DETERMINE RIGHTMOST CHARACTER TO PLOT AS '=' ON COMPLETION OR '>' OTHERWISE
if value==stop
    clear('oldpos');
    lastchar='=';
    if strcmpi(keep,'no')
        fprintf(repmat('\b',1,width+1+8));
        return;
    end
else
    oldpos=pos;
    lastchar='>';
end

% INITIALIZE OR UPDATE PROGRESS BAR
if value==start
    tstart = now();
    fprintf(1, ['[' repmat(' ',1,width-2) ']??.??.??\n']);
else
    fprintf(1, [repmat('\b',1,width+8) ...
        repmat('=',1,pos-1) lastchar repmat(' ',1,width-pos-2) ...
        ']%s\n'],...
        datestr(tstart+(now()-tstart)/(value-start)*(stop-start), 'HH.MM.SS'));
end
end


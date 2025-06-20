function draggable(h,varargin)

if length(h) > 1
    for k = 1:length(h)
        draggable(h(k),varargin{:});
    end
    return
end

% Initialization of some default arguments
user_renderer = 'zbuffer';
user_movefcn = [];
user_rotatefcn = [];
constraint = 'none';
p = [];
user_endfcn = [];       % added by SMB (see 'for k' loop below)

% At least the handle to the object must be given
Narg = nargin;
if Narg == 0
    error('Not engough input arguments');
elseif numel(h)>1
    error('Only one object at a time can be made draggable');
end

% Fetching informations about the parent axes
axh = get(h,'Parent');
if iscell(axh)
    axh = axh{1};
end
%fgh = get(axh,'Parent'); % This fails if the axes are embedded in a Panel
fgh = gcbf; % This should always work

% Assigning optional arguments
Noptarg = Narg - 1;
if Noptarg > 0
    user_movefcn = varargin{1};
    user_rotatefcn = varargin{2};
end
% =========================================================================
% Saving initial state and parameters, setting up the object callback
% =========================================================================

% Saving object's and parent figure's initial state
%getAsyncKeyState(VirtualKeyCode.VK_CONTROL);
setappdata(h,'initial_userdata',get(h,'UserData'));
setappdata(h,'initial_objbdfcn',get(h,'ButtonDownFcn'));
setappdata(h,'initial_renderer',get(fgh,'Renderer'));
setappdata(h,'initial_wbdfcn',get(fgh,'WindowButtonDownFcn'));
setappdata(h,'initial_wbufcn',get(fgh,'WindowButtonUpFcn'));
setappdata(h,'initial_wbmfcn',get(fgh,'WindowButtonMotionFcn'));

% Saving parameters
setappdata(h,'constraint_type',constraint);
setappdata(h,'constraint_parameters',p);
setappdata(h,'user_movefcn',user_movefcn);
setappdata(h,'user_rotatefcn',user_rotatefcn);
setappdata(h,'user_endfcn',user_endfcn);
setappdata(h,'user_renderer',user_renderer);

% Setting the object's ButtonDownFcn
set(h,'ButtonDownFcn',@click_object);

% =========================================================================
% FUNCTION click_object
%   Executed when the object is clicked
% =========================================================================

function click_object(obj,eventdata)
% obj here is the object to be dragged and gcf is the object's parent
% figure since the user clicked on the object
if eventdata.Button == 1
    setappdata(obj,'initial_position',get_position(obj));
    setappdata(obj,'initial_extent',compute_extent(obj));
    setappdata(obj,'initial_point',get(gca,'CurrentPoint'));
    set(gcf,'WindowButtonDownFcn',{@activate_movefcn,obj});
    set(gcf,'WindowButtonUpFcn',{@deactivate_movefcn,obj});
    activate_movefcn(gcf,eventdata,obj);
elseif eventdata.Button ~= 1
    setappdata(obj,'initial_position',get_position(obj));
    setappdata(obj,'initial_extent',compute_extent(obj));
    setappdata(obj,'initial_orientation',get_orientation(obj));
    setappdata(obj,'initial_point',get(gca,'CurrentPoint'));
    set(gcf,'WindowButtonDownFcn',{@activate_rotatefcn,obj});
    set(gcf,'WindowButtonUpFcn',{@deactivate_rotatefcn,obj});
    activate_rotatefcn(gcf,eventdata,obj);
end
% ==============================================================================
% FUNCTION activate_rotatefcn
%   Activates the WindowButtonMotionFcn for the figure
% ==============================================================================
function activate_rotatefcn(obj,eventdata,h)
set(obj,'WindowButtonMotionFcn',{@rotatefcn,h});

% ==============================================================================
% FUNCTION deactivate_rotatefcn
%   Deactivates the WindowButtonMotionFcn for the figure
% ==============================================================================
function deactivate_rotatefcn(obj,eventdata,h)
% obj here is the figure containing the object
set(obj,'WindowButtonMotionFcn',getappdata(h,'initial_wbmfcn'));
set(obj,'WindowButtonDownFcn',getappdata(h,'initial_wbdfcn'));
set(obj,'WindowButtonUpFcn',getappdata(h,'initial_wbufcn'));



% =========================================================================
% FUNCTION activate_movefcn
%   Activates the WindowButtonMotionFcn for the figure
% =========================================================================
function activate_movefcn(obj,eventdata,h)
% We were once setting up renderers here. Now we only set the movefcn
set(obj,'WindowButtonMotionFcn',{@movefcn,h});

% =========================================================================
% FUNCTION deactivate_movefcn
%   Deactivates the WindowButtonMotionFcn for the figure
% =========================================================================

function deactivate_movefcn(obj,eventdata,h)
% obj here is the figure containing the object
% Setting the original MotionFcn, DuttonDownFcn and ButtonUpFcn back
set(obj,'WindowButtonMotionFcn',getappdata(h,'initial_wbmfcn'));
set(obj,'WindowButtonDownFcn',getappdata(h,'initial_wbdfcn'));
set(obj,'WindowButtonUpFcn',getappdata(h,'initial_wbufcn'));
% Executing the user's drag end function
user_endfcn = getappdata(h,'user_endfcn');%
%getAsyncKeyState(VirtualKeyCode.VK_CONTROL);
if ~isempty(user_endfcn)
    feval(user_endfcn,h);           % added by SMB, modified by FB
end

% =========================================================================
% FUNCTION set_initial_state
%   Returns the object to its initial state
% =========================================================================

function set_initial_state(h)
initial_objbdfcn = getappdata(h,'initial_objbdfcn');
initial_userdata = getappdata(h,'initial_userdata');
set(h,'ButtonDownFcn',initial_objbdfcn);
set(h,'UserData',initial_userdata);

% =========================================================================
% FUNCTION movefcn
%   Actual code for dragging the object
% =========================================================================

function movefcn(obj,eventdata,h)
% obj here is the *figure* containing the object

% Retrieving data saved in the figure
% Reminder: "position" refers to the object position in the axes
%           "point" refers to the location of the mouse pointer
initial_point = getappdata(h,'initial_point');
p = getappdata(h,'constraint_parameters');
user_movefcn = getappdata(h,'user_movefcn');

% Getting current mouse position
current_point = get(gca,'CurrentPoint');

% Computing mouse movement (dpt is [dx dy])
cpt = current_point(1,1:3);
ipt = initial_point(1,1:3);
dpt = cpt - ipt;

% Re-computing new position with modified dpt
newpos = update_position(getappdata(h,'initial_position'),dpt);
% Setting the new position which actually moves the object
dx = set_position(h,newpos);
dx = mean(dx,1);
h.UserData.parent.set_position(h.UserData.parent.position + dx);
% Calling user-provided function handle
if ~isempty(user_movefcn)
    feval(user_movefcn,h);
end

function orient = get_orientation(obj)
props = get(obj);
%if isfield(props.UserData,'orientation')
orient = props.UserData.parent.orientation;
%else
%   error('Unable to find orientation');
%end

% =========================================================================
% FUNCTION get_position
%   Return an object's position: [x y [z / w h]] or [xdata; ydata]
% =========================================================================
function pos = get_position(obj)
props = get(obj);
if isfield(props,'Position')
    pos = props.Position;
elseif isfield(props,'Vertices')
    pos = props.Vertices';
else
    error('Unable to find position');
end

% =========================================================================
% FUNCTION update_position
%   Adds dpt to a position specification as returned by get_position
% =========================================================================
function newpos = update_position(pos,dpt)
newpos = pos;
if size(pos,1) == 1 % [x y [z / w h]]
    newpos(1:3) = newpos(1:3) + dpt;
else                % [xdata; ydata]
    newpos(1,:) = newpos(1,:) + dpt(1);
    newpos(2,:) = newpos(2,:) + dpt(2);
    newpos(3,:) = newpos(3,:) + dpt(3);
end

% =========================================================================
% FUNCTION set_position
%   Sets the position of an object obj using get_position's format
% =========================================================================
function dx = set_position(obj,pos)
if size(pos,1) == 1 % 'Position' property
    set(obj,'Position',pos);
else                % 'XData/YData' properties
    dx = pos' - obj.Vertices;
    %set(obj,'Vertices',pos');
    %
    %    set(obj,'XData',pos(1,:),'YData',pos(2,:));
end

% =========================================================================
% FUNCTION compute_extent
%   Computes an object's extent for different object types;
%   extent is [x y w h]
% =========================================================================

function extent = compute_extent(obj)
props = get(obj);
if isfield(props,'Extent')
    extent = props.Extent;
elseif isfield(props,'Position')
    extent = props.Position;
elseif isfield(props,'XData')
    minx = min(props.XData);
    miny = min(props.YData);
    w = max(props.XData) - minx;
    h = max(props.YData) - miny;
    extent = [minx miny w h];
else
    error('Unable to compute extent');
end

% =========================================================================
% FUNCTION is_inside_range
%   Checks if a rectangular object is entirely inside a rectangular range
% =========================================================================

function inrange = is_inside_range(extent,range)
% extent is in the [x y w h] format
% range is in the [xmin xmax ymin ymax] format
% inrange is a 4x1 vector of boolean values corresponding to range limits
inrange = [extent(1) >= range(1) ...
    extent(1) + extent(3) <= range(2) ...
    extent(2) >= range(3) ...
    extent(2) + extent(4) <= range(4)];

% ==============================================================================
% FUNCTION rotatefcn
%   Actual code for dragging the object
% ==============================================================================

function rotatefcn(obj,eventdata,h)
% obj here is the figure containing the object

% Getting current point
current_point = get(gca,'CurrentPoint');
initial_point = getappdata(h,'initial_point');
user_rotatefcn = getappdata(h,'user_rotatefcn');



fi0 = getappdata(h,'initial_orientation');
speed = 1;
ori_act = h.UserData.parent.orientation;
switch gca().UserData
    case 'x'
        cpt = current_point(1,2:3);
        ipt = initial_point(1,2:3);
        dpt = cpt - ipt;
        setappdata(h,'old_dpt',dpt);
        dfi = sum(dpt)*speed  -  atan2(ori_act(3),ori_act(2)) + atan2(fi0(3),fi0(2));
        %origin = [0,mean(h.Vertices(:,2:3),1)];
        %rotate(h,[1 0 0], dfi*180/pi,  origin);
        M = [1, 0, 0;
            0, cos(dfi), -sin(dfi);
            0, sin(dfi), cos(dfi)];
    case 'y'
        cpt = [-current_point(1,1),current_point(1,3)];
        ipt = [-initial_point(1,1),initial_point(1,3)];
        dpt = cpt - ipt;
        setappdata(h,'old_dpt',dpt);
        dfi = sum(dpt)*speed  -  atan2(ori_act(3),ori_act(1)) + atan2(fi0(3),fi0(1));
        %origin = [mean(h.Vertices(:,1)),0,mean(h.Vertices(:,3))];
        %rotate(h,[0 1 0], dfi*180/pi,  origin);

        M = [cos(dfi),0, sin(dfi);
            0 , 1 , 0
            -sin(dfi),0, cos(dfi)];
    case 'z'
        cpt = current_point(1,1:2);
        ipt = initial_point(1,1:2);
        dpt = cpt - ipt;
        setappdata(h,'old_dpt',dpt);
        dfi = sum(dpt)*speed  -  atan2(ori_act(2),ori_act(1)) + atan2(fi0(2),fi0(1));
        %origin = [mean(h.Vertices(:,1:2),1),0];

        %rotate(h,[0 0 1], dfi*180/pi,  origin);
        M = [cos(dfi), -sin(dfi),0;
            sin(dfi), cos(dfi),0;
            0,0,1];
    case '3D'
        print('Rotation is not allowed in 3D view');
        M = eye(3);
end
h.UserData.parent.set_orientation((M* ori_act')');
if ~isempty(user_rotatefcn)
    feval(user_rotatefcn,h);
end

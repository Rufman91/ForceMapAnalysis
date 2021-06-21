%
% PURPOSE:
%
%   Show how to add drop support from file explorer to some matlab axis
%
% SYNTAX:
%
%   [] = DropListenerDemo();
%
% USAGE:
%
%   Simply drop files from file explorer into displayed axis.
%

%%
function [] = DropListenerDemo(String)
%[
    % Create a figure with some axis inside
    
    if nargin < 1
        String = '';
    end
    
    fig = figure(666); clf;
    List = uicontrol('Parent', fig,...
        'Style','listbox',...
        'Max',10000,'Min',1,...
        'Units','normalized',...
        'Position',[0.1 0.1 0.8 0.8]);
    
    set(List,'string',String)
    
    % Get back the java component associated to the axis
    % NB1: See §3.7.2 of Undocumented Secrets of Matlab Java Programming
    % NB2: or use findjobj, or javaObjectEDT for drop support onto other component types
    jFrame = get(handle(fig), 'JavaFrame');
    jAxis = jFrame.getAxisComponent();
    
    % Add listener for drop operations
    DropListener(jAxis, ... % The component to be observed
                 'DropFcn', @(s, e)onDrop(fig, s, e)); % Function to call on drop operation    
%]
end
function [] = onDrop(fig, listener, evtArg) %#ok<INUSL>
%[
    % Get back the dropped data
    data = evtArg.GetTransferableData();
    
    % Is it transferable as a list of files
    if (data.IsTransferableAsFileList)       
        
        % Do whatever you need with this list of files
                
        String = data.TransferAsFileList;
        DropListenerDemo(String)
        % Indicate to the source that drop has completed 
        evtArg.DropComplete(true);
        
    elseif (data.IsTransferableAsString)
        
        % Not interested
        evtArg.DropComplete(false);
        
    else
        
        % Not interested
        evtArg.DropComplete(false);
        
    end
%]
end
function update_suggestion_list(source,event,MatchCase,ListBox,ACDictionary)

if isequal(event.Key,'shift')
    return
end

disp('Callback triggered!')

if ~isequal(event.Key,'return')
    import java.awt.Robot;
    import java.awt.event.KeyEvent;
    robot=Robot;
    robot.keyPress(KeyEvent.VK_ENTER);
    pause(0.01)
    robot.keyRelease(KeyEvent.VK_ENTER);
end

EditText = get(source,'String');

MatchedDic = match_string_to_dictionary(EditText,ACDictionary,MatchCase);

ListBox.String = MatchedDic;
ListBox.Value = 1;

end
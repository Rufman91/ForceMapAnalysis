function autofill_edit_field_from_list(source,event,EditBox)

if isempty(source.String)
    return
end

ChosenString = source.String{source.Value};

EditBox.String = ChosenString;

uicontrol(EditBox)

end
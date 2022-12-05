function update_suggestion_list(source,event,ListBox)

disp(source.String)

set(ListBox,'Value',rem(length(split(get(source,'String'),' ')),length(ListBox.String)));

split(source.String,' ')

end
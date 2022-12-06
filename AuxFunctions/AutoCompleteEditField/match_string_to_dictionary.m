function MatchedDic = match_string_to_dictionary(String,Dic,MatchCase)

if ~MatchCase
    PrepString = lower(String);
    PrepDic = lower(Dic);
else
    PrepString = String;
    PrepDic = Dic;
end

Matching = contains(PrepDic,PrepString);

MatchedDic = Dic(Matching);

end
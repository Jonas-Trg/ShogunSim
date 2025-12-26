function WordCaption(ActXWord,type,text,enter)
% WordCaption
% 1.0.0.0
if ishandle(ActXWord)
    if any(strcmp(type,{'Table','Figure'}))
        ActXWord.Selection.InsertCaption(type,[': ' text],'',1)
        if enter(end)
            WordText(ActXWord,'',[],[0 1]);
        end
    end
end
return %WordCaption

function WordIndent(actx_word,cm_left,cm_firstline,style)
% WordIndent
% 1.0.0.0
if ishandle(actx_word)
    if nargin >= 4 && ~isempty(style)
        if strcmp(style,'Underline')
            actx_word.Selection.Font.Underline=1;
        elseif strcmp(style,'Bold')
            actx_word.Selection.Font.Bold=1;
        elseif strcmp(style,'Italic')
            actx_word.Selection.Font.Italic=1;
        else
            actx_word.Selection.Style = style;
        end
    end
    actx_word.Selection.ParagraphFormat.LeftIndent=cm_left*28.35;
    actx_word.Selection.ParagraphFormat.FirstLineIndent=cm_firstline*28.35;
end
return %WordIndent

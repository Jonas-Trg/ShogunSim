function [actx_word,word_handle]=WordStart(word_file_p)
% WordStart
% 1.0.0.0
% Start an ActiveX session with Word:
actx_word = actxserver('Word.Application');
actx_word.Visible = true;
trace(actx_word.Visible);
if ~exist(word_file_p,'file')
    % Create new document:
    word_handle = invoke(actx_word.Documents,'Add');
else
    % Open existing document:
    word_handle = invoke(actx_word.Documents,'Open',word_file_p);
end

return % WordStart

function WordGoToBkMrk(ActXWord,BkMrk)
% WordGoToBkMrk
% 1.0.0.0
if ishandle(ActXWord)
    ActXWord.ActiveDocument.Bookmarks.Item(BkMrk).Select
else
    fprintf([BkMrk ': '])
end
return %WordGoToBkMrk

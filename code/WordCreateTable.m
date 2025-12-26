function WordCreateTable(actx_word,cellData,style,enters)
% WordCreateTable
% 1.0.0.0
if ishandle(actx_word)
    tSize=size(cellData);
    %Add a table which auto fits cell's size to contents
    if(enters(1))
        actx_word.Selection.TypeParagraph; %enter
    end
    %create the table
    %Add = handle Add(handle, handle, int32, int32, Variant(Optional))
    actx_word.ActiveDocument.Tables.Add(actx_word.Selection.Range,tSize(1),tSize(2),1,1);
    %Hard-coded optionals
    %first 1 same as DefaultTableBehavior:=wdWord9TableBehavior
    %last  1 same as AutoFitBehavior:= wdAutoFitContent

    if size(style,2)~=tSize(2)
        if size(style,2) == 1
            style(1,2:tSize(2))=style(1,1);
        else
            style{1,tSize(2)}=[];
        end
    end
    for r=1:tSize(1)
        for c=1:tSize(2)
            if r==1
                WordText(actx_word,cellData{r,c},{'Bold',style{c}},[],[],1);
                actx_word.Selection.Rows.HeadingFormat=-1;
            else
                WordText(actx_word,cellData{r,c},style{c});
            end
            actx_word.Selection.MoveRight;
        end
        actx_word.Selection.MoveRight;
    end
else
    tSize=0;
    for ii=1:size(cellData,1)
        for jj=1:size(cellData,2)
            tSize=max(tSize,length(cellData{ii,jj})+1);
        end
    end
    tab=ceil(tSize/4)*4;
    
    for ii=1:size(cellData,1)
        for jj=1:size(cellData,2)
            str=strrep(cellData{ii,jj},'\_','_');
            fprintf(pad(str,tab))
        end
        fprintf(newline)
    end
end
return % WordCreateTable



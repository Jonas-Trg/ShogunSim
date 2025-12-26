function WordText(actx_word,txt,style,enters,colors,keepWithNext)
% WordText
% 1.0.0.0
%
%   input
%       actx_word - com object to MS Word
%       txt - text to be entered
%       style (optional) - reference to build in Styles in Word doc
%       enters (optional) - 1x2 of number of "enters" before and after text
%       colors (optional) - text colour
%       keepWithNext (optional) - no page break btw this and next paragr
%   Description
%       txt can include formatting characters "^" for superscript and "_"
%       for subscript. Example: 'Acc: 1.0 m/s^2' - '2' will be superscript
%       'F_{Traction} = 10 kN' - 'Traction' will be subscript
%
if ishandle(actx_word)
    if (nargin >= 4) && ~isempty(enters) && (enters(1))
        for ii=1:enters(1)
            actx_word.Selection.TypeParagraph; %enter
        end
    end
    if nargin >= 3 && ~isempty(style)
        if ~iscell(style)
            style={style};
        end
        for ii=1:size(style,2)
            if strcmpi(style{ii},'Underline')
                actx_word.Selection.Font.Underline=1;
            elseif strcmpi(style{ii},'Bold')
                actx_word.Selection.Font.Bold=1;
            elseif strcmpi(style{ii},'Italic')
                actx_word.Selection.Font.Italic=1;
            elseif strcmpi(style{ii},'Left')
                actx_word.Selection.ParagraphFormat.Alignment=0;
            elseif strcmpi(style{ii},'Centered')
                actx_word.Selection.ParagraphFormat.Alignment=1;
            elseif strcmpi(style{ii},'Right')
                actx_word.Selection.ParagraphFormat.Alignment=2;
            elseif ~isempty(style{ii})
                actx_word.Selection.Style = style{ii};
            end
        end
    end
    if(nargin >= 5) && ~isempty(colors) %check to see if color_p is defined
        actx_word.Selection.Font.ColorIndex=colors;
    end
    b=find(txt=='_' | txt=='^',1,'first');
    while ~isempty(b)
        if b>1 && txt(b-1)~='\'
            actx_word.Selection.TypeText(txt(1:b-1));
            if txt(b)=='_'
                actx_word.Selection.Font.Subscript=1;
            elseif txt(b)=='^'
                actx_word.Selection.Font.Superscript=1;
            end
            if txt(b+1)=='{'
                c=find(txt(b+2:end)=='}',1,'first')+b;
                actx_word.Selection.TypeText(txt(b+2:c));
                txt=txt(c+2:end);
            else
                actx_word.Selection.TypeText(txt(b+1));
                txt=txt(b+2:end);
            end
            actx_word.Selection.Font.Superscript=0;
            actx_word.Selection.Font.Subscript=0;
            b=find(txt=='_' | txt=='^',1,'first');
        else
            actx_word.Selection.TypeText(txt([1:b-2 b]));
            txt=txt(b+1:end);
            b=find(txt=='_' | txt=='^',1,'first');
        end
    end
    if isempty(txt)
        txt='';
    end
    actx_word.Selection.TypeText(txt);
    if(nargin >= 5) && ~isempty(colors) %check to see if color_p is defined
        actx_word.Selection.Font.ColorIndex='wdAuto';%set back to default color
    end
    if(nargin == 6) && ~isempty(keepWithNext)
        actx_word.Selection.ParagraphFormat.keepWithNext=-keepWithNext;
    end
    if nargin >= 4 && size(enters,2)>1
        for k=1:enters(2)
            actx_word.Selection.TypeParagraph; %enter
        end
    end
    if nargin >= 3 && ~isempty(style)
        if strcmp(style,'Underline')
            actx_word.Selection.Font.Underline=0;
        elseif strcmp(style,'Bold')
            actx_word.Selection.Font.Bold=0;
        elseif strcmp(style,'Italic')
            actx_word.Selection.Font.Italic=0;
        end
    end
elseif size(txt,2)>0
    txt=strrep(txt,'\_','_');
    txt=strrep(txt,'\','\\');
    if (nargin >= 4) && ~isempty(enters) && (enters(1))
        fprintf(newline)
    end
    b=find(txt==char(9),1,'first');
    fprintf(pad(txt(1:b),50));
    fprintf(txt(b+1:end));
    if (nargin >= 4) && ~isempty(enters) && (enters(2))
        fprintf(newline)
    end
end
return % WordText
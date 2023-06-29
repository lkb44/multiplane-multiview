function mosaic4D(varargin)
% mosaic creates subplots with the data in the different fields 
% 

if isa(varargin{end},'function_handle')
    feval(varargin{end});
    return;
else
    
    %
    % get the different matrices
    %
    udata.nargin = nargin;
    if nargin>=1
        udata.mats = varargin{1};
        disp(['Display data: ' num2str(size(udata.mats)) ] )
        if (nargin>=2)
            udata.masks = varargin{2};
            disp(['Use mask of size: ' num2str(size(udata.masks)) ] )
        end
    end
    
    figure;
    
    slider = uicontrol('Style', 'slider',...
       'Max',size(udata.mats,3),'Min',1,...
       'Units', 'normalized', ...
       'Position', [.25 .005 .4 .04],...
       'SliderStep', [1/(size(udata.mats,3)-1) 1/(size(udata.mats,3)-1)], ...
       'Value', 1, ...
       'Callback', {@mosaic4D,@move_slider});
    text_box = uicontrol('Style', 'text',...
       'Units', 'normalized', ...
       'Position', [.675 .006 .04 .04],...
       'String', '1');
   
    udata.text_box = text_box;
    set(gcf,'UserData',udata);
        showimg(1, udata);
end

end


function move_slider
udata = get(gcf,'UserData');
val = round(get(gcbo,'Value'));
set(gcbo,'Value',val);
set(udata.text_box,'String',num2str(val));
save_xlim = xlim;
save_ylim = ylim;
showimg(val,udata);
xlim(save_xlim);
ylim(save_ylim);
end

%where the images are displayed 
function showimg(idx,udata)
    mats = udata.mats;
    nargin = udata.nargin;
    if nargin>=2
        masks = udata.masks;
    end
    planes = size(mats,4);
    for i= 1:planes
      names{i} = i;
    end
    dim = ceil(sqrt(planes));
    
 %%%%%%---------HOW SHOULD 2D IMAGE BE TRANSFORMED-------%
 transfim = @(x) x; %transpose(x);
    for matNum = 1:planes ; 
        subplot(ceil(planes/dim),dim , matNum) 
        imdisp(transfim(mats(:,:,idx,matNum)));
%         imdisp((mats(:,:,idx,matNum)));

        if (nargin>=2)
        hold on;
        s = sumall( masks(:,:,idx).^2 );
        if s>0
            contour(transfim(masks(:,:,idx)),'Color','yellow');
%             contour((masks(:,:,idx)),'Color','yellow');

%             disp(['Cant display mask, null data'])
        end
        hold off
        end
        title(['slice' num2str(names{matNum}) ' at timepoint ' num2str(idx) ])
    end
end

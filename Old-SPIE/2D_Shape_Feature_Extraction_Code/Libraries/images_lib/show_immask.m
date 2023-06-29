function show_immask(im,mask)
       
 
    im = rescale(im);
    
    if ~isempty(mask)
        if ~isempty(find(mask>0,1))
           colors = [0 1 0; 1 1 0; 1 2/3 0; 0 0 1; 1 0 0; 1 0 1; 0 1 1; 1 1 1;0 0 0]; %[green yellow orange blue red magenta cyan white black]
           if length(find(unique(mask)>0)) > length(colors)
                error(['too many label colors! (max allowed = ', int2str(length(colors)),', currently)']);
           elseif max(unique(mask)) > length(colors)
               error(['Color value of a label is greater than current number of color bins!']);
           end

            %format the order in such a way that the biggest masks are displayed first, so that you can see all colors:)
            uvec = mask(mask>0);
            if ~isempty(uvec)
                u = unique(uvec);
                N = histcounts(uvec,length(u));
                [~,ind] = sort(N,'descend'); %*******descend, ascend (whichever helps see all colors!)********
                for i = 1:length(ind)
                    m = (mask==u(ind(i)));
                    [dx,dy] = gradient(m);
                    edgei=sqrt(dx.^2+dy.^2);
                    edgei=edgei~=0;
    %                 se=[0 1 0; 1 1 1; 0 1 0];
    %                 se=se~=0;
    %                 edgei=imerode(edgei,se);
                    im = rgbmaskrgb(im,edgei,colors(u(ind(i)),:));
                end
            end
            imagesc(im);colormap gray; axis image;axis off;hold on;
        else
            imagesc(im);axis image;axis off;colormap gray;
        end
    end
      
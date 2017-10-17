function [X,Y,u0,M,N] = Interactive
    close all
    f = figure(2);
    str = 'Size: 3';
    start = uicontrol('Style','pushbutton',...
                        'Position',[40,300,40,20],...
                        'String','Start',...
                        'Callback',@startcall);

    clr = uicontrol('Style','pushbutton',...
                    'Position',[40,200,40,20],...
                    'String','Clear',...
                    'Callback',@clearcall);
                
    slr_anno = annotation('textbox',[.05,.27,.1,.05],...
                          'String',str);
                    
    slider = uicontrol('Style','slider',...
                        'Min',1,'Max',9,'Value',3,...
                        'SliderStep',[1/4,1/2],...
                        'Position',[20,150,80,20],...
                        'Callback',{@barcall,slr_anno});

                    
                    
    im1 = imread('cameraman.jpg');
    im1 = rgb2gray(im1);
    u0 = im2double(im1);
    N = size(u0,2); % # of columns, x direction;
    M = size(u0,1); % # of rows, y direction;
    imshow(u0)
    xlim([0,N]);
    ylim([0,M]);

    Flag = 1;
    X = [];
    Y = [];
    w =3;
    cursor_handle = rectangle('Position',[100,100,w,w],...
                                'FaceColor','r',...
                                'EdgeColor','r',...
                                'visible','off');
    function startcall(gcbo,eventdata)
        cursor_handle = [];
        uiresume(gcf)
    end


    function clearcall(gcbo,eventdata)
        X = [];
        Y = [];
        f = imshow(im1);
    end

    function barcall(gcbo,eventdata,text_handle)
        
        w = get(gcbo,'Value');
        str = ['Size: ' num2str(w)];
        set(text_handle,'String',str);
        
    end



    function mouseclick(gcbo,eventdata)
     try
        Flag = 0;
        cp = get(gca, 'Currentpoint');
        x = round(cp(1,1));
        y = round(cp(1,2));
        pos_x = x-(w-1)/2;
        pos_y = y-(w-1)/2;
        set(cursor_handle,'Position',...
                           [pos_x,pos_y,w,w],...
                           'visible','on');
        rectangle('Position',[pos_x,pos_y,w,w],...
                   'FaceColor','r','EdgeColor','r')
        x = pos_x:(pos_x+w-1);
        y = pos_y:(pos_y+w-1);
        [x,y] = meshgrid(x,y);
        x = x(:);
        y = y(:);
        X = [X;x];
        Y = [Y;y];
     end 
               
        
    end

    function mousemove(gcbo,eventdata)
        if Flag == 0
           try
            cp = get(gca, 'Currentpoint');
            x = round(cp(1,1));
            y = round(cp(1,2));
            pos_x = x-(w-1)/2;
            pos_y = y-(w-1)/2;
            set(cursor_handle,'Position',...
                           [pos_x,pos_y,w,w],...
                            'visible','on');
            rectangle('Position',[pos_x,pos_y,w,w],...
                    'FaceColor','r','EdgeColor','r')
            x = pos_x:(pos_x+w-1);
            y = pos_y:(pos_y+w-1);
            [x,y] = meshgrid(x,y);
            x = x(:);
            y = y(:);
            X = [X;x];
            Y = [Y;y];
           end
        end
    end

    function mouserealse(gcbo,eventdata)
        hold on
        Flag = 1;
        try
            set(cursor_handle,'visible','off');
        end
    
    end

    set(f,'WindowButtonDownFcn', @mouseclick);
    set(f,'WindowButtonMotionFcn',@mousemove);
    set(f,'WindowButtonUpFcn',@mouserealse);
    uiwait(gcf)
    

end
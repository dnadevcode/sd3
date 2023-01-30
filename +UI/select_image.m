function varargout = select_image(listing,output,ii,nameMolAll,idxMovAll,dim1,dim2)
    h=figure('CloseRequestFcn',@my_closereq,'units','normalized','outerposition',[0 0 1 1])
%     set(gca,'XTick',[])
%                     set(gca,'XColor','none')

    iiSt = 1;
    import UI.plot_mol;
    
    tiledlayout(dim1,2,'TileSpacing','none','Padding','tight');

    try
        h1 = [];
        for idx = 1:dim1
            for jdx = 1:dim2                
                file = fullfile(listing(iiSt+ii-1).folder,listing(iiSt+ii-1).name);
                NN = length(imfinfo(file));
                img = cell(1,NN);
                for iiD=1:NN % Bars, Dots, DotsM, Bitmask, BarsM
                    img{iiD} = imread(file,iiD);
                end
                % img1

%                 spltName = strsplit(listing(iiSt+ii-1).name,'_mol_');
%                 nameMov = spltName{1};
%                 spltName2 = strsplit(spltName{2},'.');
%                 nameMol = str2num(spltName2{1});
%                 names = cellfun(@(x) x.name,output,'un',false);
%                 idxMov = find(cellfun(@(x) ~isempty(strfind(x,nameMov)),names));
                nameMol = nameMolAll(iiSt+ii-1);
                idxMov = idxMovAll(iiSt+ii-1);

                nexttile;
                h1(iiSt) = imagesc(img{1});
                vOff = output{idxMov}.boundaries{nameMol}(1);
                hOff = output{idxMov}.boundaries{nameMol}(3);
                str = sprintf('Mol. %i',nameMol);

                text(hOff,-5+vOff,str,'Color','white');%
                set(gca,'XColor','none')
                set(gca,'YColor','none')

                nexttile
                h1(iiSt) = imagesc(img{2});
                hold on;
                plot_mol(output,idxMov,nameMol);
                plot(-output{idxMov}.pos{nameMol}(2)+output{idxMov}.trueedge{nameMol}(:,2)+1,-output{idxMov}.pos{nameMol}(1)+output{idxMov}.trueedge{nameMol}(:,1)+1,'black.')
                set(gca,'XColor','none');
                set(gca,'YColor','none');
                set(h1(iiSt), 'buttondownfcn', {@loads_of_stuff,iiSt+ii-1});
                iiSt = iiSt+1;
            end
        end
    catch
    end
    uiwait()
    
    
    

    function loads_of_stuff(src,eventdata,x)
        if get(src,'UserData')
            set(src,'UserData',0)
%             title('');
            text(1,1,'Removed','Color','white','BackgroundColor','blue');

        else
            set(src,'UserData',1)
%             title('Selected');
            text(1,1,'Selected','Color','white','BackgroundColor','blue');
        end
%         fprintf('%s\n',num2str(x));
%         C = get(h, 'UserData')
    
    end
%     
function my_closereq(src,callbackdata)
% Close request function 
% to display a question dialog box 
    try
    varargout{1} = find(cellfun(@(x) ~isempty(x),get(h1,'Userdata')));
    catch
    end
    delete(h)
%     uiresume() 

end

end

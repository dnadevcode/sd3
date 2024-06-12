function [goodMols] = pairwise_distances_molecule(trueedge,minSeparation)
    % Finds molecules within minSeparation from each other

goodMols = true(1,length(trueedge));
% tic
%     for i=1:length(trueedge)-1
%         v1x=unique(trueedge{i}(:,1));
% %         v1y=unique(trueedge{i}(:,2));
% 
%         for j=i+1:length(trueedge)
% %             trueedge{i}        
%             v2x=unique(trueedge{j}(:,1));
% %             v2y=unique(trueedge{j}(:,2));
%                 
% %             [k,dist] = dsearchn(trueedge{i},trueedge{j});
%             [Idx,md] = knnsearch(trueedge{i}, trueedge{j});
%             if min(md) <= minSeparation % molecules to remove
%                 goodMols([i j]) = false;
%             end
%         end
%     end


    for i=1:length(trueedge)-1
        for j=i+1:length(trueedge)
%              trueedge{i}
                
%             [k,dist] = dsearchn(trueedge{i},trueedge{j});
            [Idx,md] = knnsearch(trueedge{i}, trueedge{j});
            if min(md) <= minSeparation % molecules to remove
                goodMols([i j]) = false;
            end
        end
    end
% toc
% 
%%
% i=6
% j=3
% figure
% plot(trueedge{i}(:,1),trueedge{i}(:,2),'ko')
% hold on
% plot(trueedge{j}(:,1),trueedge{j}(:,2),'*g')
%%
end


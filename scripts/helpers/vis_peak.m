function[fig] = vis_peak(stack,mat,peak,pad_xy,pad_z,main_title)

num_channels = size(stack,4);
num_cycles = size(stack,5);

xleft = max(peak(:,1)-pad_xy,1); xright = min(peak(:,1)+pad_xy,size(stack,1));
yleft = max(peak(:,2)-pad_xy,1); yright = min(peak(:,2)+pad_xy,size(stack,2));
zleft = max(peak(:,3)-pad_z,1); zright = min(peak(:,3)+pad_z,size(stack,3));


for cycle=1:num_cycles
    for channel=1:num_channels
        zoom_stack(:,:,channel,cycle) = squeeze(max(stack(xleft:xright,yleft:yright,zleft:zright,channel,cycle),[],3));   
    end
end

fig = figure('visible','off');
p = tight_subplot(num_channels+1,num_cycles,[0.005 0.005],[0.005 0.005],[0.005 0.005]);

for cycle=1:num_cycles
    
    cy_max = max(reshape(zoom_stack(:,:,:,cycle),[],1));
    cm = jet(floor(cy_max));
    
    for channel=1:num_channels+1
        
        %axes(p((channel-1)*num_cycles+cycle)); 
        axes(p((channel-1)*num_cycles+cycle)); 
        
        if channel<=4

        imshow(capImage(zoom_stack(:,:,channel,cycle),cy_max,'abs'),cm); hold on;
        %imshow(capImage(zoom_stack(:,:,channel,cycle),cy_max,'abs'),jet); hold on;
        rectangle('Position',[pad_xy-1,pad_xy-1,5,5],'EdgeColor','magenta');
        title(sprintf('%.02f',mat(channel,cycle).^2./sum(mat(:,cycle).^2,1)))
        
        else
        
        imshow(capImage(max(zoom_stack(:,:,:,cycle),[],3),cy_max,'abs'),cm); hold on;
        %imshow(capImage(zoom_stack(:,:,channel,cycle),cy_max,'abs'),jet); hold on;
        rectangle('Position',[pad_xy-1,pad_xy-1,5,5],'EdgeColor','magenta');
        
        end
        
    end
    
end

st = suptitle(main_title)
set(st,'FontSize',10,'FontWeight','normal')

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 num_cycles num_channels+2.5];

end
%% for making RGB vis

%{
    cap = prctile(reshape(zoom_stack(:,:,:,cycle),[],1),99);
        
    for channel=1:num_channels
        axes(p((channel-1)*num_cycles+cycle)); 
        rgb_image = zeros(size(zoom_stack,1),size(zoom_stack,2),3);
        if channel == 1
            rgb_image(:,:,1) = zoom_stack(:,:,channel,cycle)./max(reshape(zoom_stack(:,:,:,cycle),[],1)).*1;
        elseif channel == 2
            rgb_image(:,:,1) = zoom_stack(:,:,channel,cycle)./max(reshape(zoom_stack(:,:,:,cycle),[],1)).*1;
            rgb_image(:,:,3) = zoom_stack(:,:,channel,cycle)./max(reshape(zoom_stack(:,:,:,cycle),[],1)).*1;
        elseif channel == 3
            rgb_image(:,:,2) = zoom_stack(:,:,channel,cycle)./max(reshape(zoom_stack(:,:,:,cycle),[],1)).*1;
        elseif channel == 4
            rgb_image(:,:,3) = zoom_stack(:,:,channel,cycle)./max(reshape(zoom_stack(:,:,:,cycle),[],1)).*1;
        end
        imshow(rgb_image,[]); hold on;
        %imshow(region_stack(:,:,channel,cycle),[0 cap]); hold on;
    end

%}
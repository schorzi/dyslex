%%
close all; clear variables;
% set parameters
white_val             = 255;
frame_size            = 100;
max_image_size_factor = 0.8;
num_of_sizes          = 3;
num_of_angles         = 6;

% organize data
num_of_images         = num_of_sizes*num_of_angles;
object_max_size       = frame_size*max_image_size_factor;

raw_image      = imread('thumbup.png');
shape_max_size = imresize(raw_image, [object_max_size object_max_size]);
pad_base       = ceil(0.5*frame_size*(1-max_image_size_factor));
base_image     = padarray(shape_max_size,[pad_base pad_base],white_val);

sizes_vec         = 2*round([frame_size:-(0.75*frame_size/num_of_sizes):frame_size/4]/2);
angles_vec        = 0:360/num_of_angles:360-(360/num_of_angles);
rand_sizes_order  = randperm(num_of_sizes);
rand_angles_order = randperm(num_of_angles);
data              = zeros(frame_size,frame_size,num_of_images);

%  create diverse images
for i=1:num_of_sizes   
    shape_resize = imresize(base_image, [sizes_vec(rand_sizes_order(i)) sizes_vec(rand_sizes_order(i))]);
    for j=1:num_of_angles
        shape_rotated  = imrotate(shape_resize,angles_vec(rand_angles_order(j)),'nearest','crop');
        Mrot           = ~imrotate(true(size(shape_resize)),angles_vec(rand_angles_order(j)),'nearest','crop');
        shape_rotated(Mrot&~imclearborder(Mrot)) = white_val;
        new_size      = size(shape_rotated,1);
        pad_num       = ceil((frame_size-new_size)/2);
        data(:,:,(i-1)*num_of_angles+j)   = padarray(shape_rotated,[pad_num pad_num],white_val);
    end
    
end
% show image montage 
 data_montage(:,:,1,:) = data; 
 montage(data_montage); title('Data Images');
%% run diffusion-map algorithm
mX = reshape(data,[],num_of_sizes*num_of_angles);
mW  = squareform( pdist(mX') );
eps = .5 * median(mW(:).^2);
mK  = exp(- mW.^2 / eps);
mA  = bsxfun(@rdivide, mK, sum(mK, 2));

[mU, mL, mV] = eig(mA);

mP1 = mU(:, 2);
mP2 = mU(:, 2:3);
mP3 = mU(:, 2:4);

%% plot results
figure; plot(mP1); title ('1st eig-vec');
figure; scatter(mP2(:,1), mP2(:,2), 20, 'Fill');  title('first 2 eig-vecs');
figure; scatter3(mP3(:,1), mP3(:,2), mP3(:,3), 20, 'Fill'); title 'first 3 eig-vecs');

%%
figure; 
minZYdist = sqrt(min(pdist(mP3(:,[2 3]))));
coordinate_resize_factor = 5*frame_size/minZYdist;
for i=1:num_of_images
    img_resize = imresize(data(:,:,i), 1);
    
    x = mP3(i,1)*coordinate_resize_factor;
    y = mP3(i,2)*coordinate_resize_factor;
    z = mP3(i,3)*coordinate_resize_factor;
    
    xImage = [x x; x x];
    yImage = [y+(frame_size/2) y+(frame_size/2); y-(frame_size/2) y-(frame_size/2)];
    zImage = [z-(frame_size/2) z+(frame_size/2); z-(frame_size/2) z+(frame_size/2)];

    surf(xImage,yImage,zImage,'CData',img_resize,...
     'FaceColor','texturemap');
    hold on
end
title('first 3 eig-vecs');
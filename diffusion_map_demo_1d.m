%%
close all; clear variables;

white_val             = 255;
frame_size            = 100;
max_image_size_factor = 0.8;
num_of_sizes          = 20

object_max_size       = frame_size*max_image_size_factor;

raw_image      = imread('thumbup.png');
shape_max_size = imresize(raw_image, [object_max_size object_max_size]);
pad_base       = ceil(0.5*frame_size*(1-max_image_size_factor));
base_image     = padarray(shape_max_size,[pad_base pad_base],white_val);

sizes_vec  = frame_size:-(0.8*frame_size/num_of_sizes):frame_size/5;
rand_order = randperm(num_of_sizes);
data       = zeros(frame_size,frame_size,num_of_sizes);

for i=1:num_of_sizes
    
    shape_resize = imresize(base_image, [sizes_vec(rand_order(i)) sizes_vec(rand_order(i))]);
    new_size     = size(shape_resize,1);
    pad_num      = ceil((frame_size-new_size)/2);
    data(:,:,i)      = padarray(shape_resize,[pad_num pad_num],white_val);
end

 data_montage(:,:,1,:) = data; 
 montage(data_montage);
%%
mX = reshape(data,[],num_of_sizes);
mW  = squareform( pdist(mX') );
eps = .5 * median(mW(:).^2);
mK  = exp(- mW.^2 / eps);
mA  = bsxfun(@rdivide, mK, sum(mK, 2));

[mU, mL, mV] = eig(mA);

mP1 = mU(:, 2);
mP2 = mU(:, 2:3);
mP3 = mU(:, 2:4);

figure; plot(mP1);
figure; scatter(mP2(:,1), mP2(:,2), 100, mW(1,:)', 'Fill');   colorbar;
figure; scatter3(mP3(:,1), mP3(:,2), mP3(:,3), 100, 'Fill'); colorbar;

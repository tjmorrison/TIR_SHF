%% Make pcolor

%Camera specs
cam_height=1.67; %meters
cam_angle_horizon=29.23; %degrees
angle_y=21.7; %degrees
angle_x=17.5; %degrees
npixel_x=128;
npixel_y=160;

[pixel_x,pixel_y]=camera_footprint_height_angle(cam_height,cam_angle_horizon,...
    angle_y,angle_x,npixel_y,npixel_x);

figure();
pcolor(pixel_x, pixel_y,flip(squeeze(mean(q_air(1:128,1:160,:),3)')));
shading interp; 
colorbar;
xlabel('x (m)');
ylabel('y (m)')
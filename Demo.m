%% Demo for Frequency Boundary Energy Model 
% This code contain 3 demos.  
demo =2; % 1,2,3
switch demo
    case 1
        Domain  = 0; % using FFT to perform FBEM
        max_its = 700; %
        length  = 0.0015; % for length term
        regul   = 0.1; % for level set function regularization term
        thresh  = 2; % terminal condition
        for i =1:5     
            if i<4
                img = imread(['demo1_1.tif']);
            else
                img = imread(['demo1_',int2str(i-2),'.tif']);
            end
            load(['demo1_',int2str(i),'.mat'])
            [seg,phi,its] = fbe_acm(double(img),init_mask,max_its,length,regul,...
            sigma1,sigma2,Domain,thresh);
        end
    case 2
         for i =1:1   
            img = imread(['demo2_',int2str(i),'.jpg']);
            load(['demo2_',int2str(i),'.mat'])
            max_its = 30;
            [seg,phi,its] = fbe_acm(double(img),init_mask,max_its,length,regularization,...
            sigma1,sigma2,domain,thresh);%,'r',false);
        end
    case 3
         for i =1:2     
            img = imread(['demo3_',int2str(i),'.jpg']);
            load(['demo3_',int2str(i),'.mat'])
            [seg,phi,its] = fbe_acm(double(img),init_mask,max_its,length,regularization,...
            sigma1,sigma2,domain,thresh,'r',false);
        end        
end
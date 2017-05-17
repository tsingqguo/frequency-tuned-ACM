
function [seg,phi,its] = fbe_acm(img,init_mask,max_its,length,regularization,...
    sigma1,sigma2,Domain,thresh,color,display)
    
    if(~exist('max_its','var')) 
        max_its = 200; 
    end
    %-- default value for parameter length is 1
    if(~exist('length','var')) 
        length = 1; 
    end
    %-- default value for parameter penalizing is 1
    if(~exist('regularization','var')) 
        regularization = .1; 
    end    
    %-- default value for parameter sigma1
    if(~exist('sigma1','var')) 
        sigma1 = 9;
    end
    %-- default value for parameter sigma2
    if(~exist('sigma2','var')) 
        sigma2 = 3;
    end
    % Domain =0 ,using FFT, otherwise, using spatial filter
    if(~exist('Domain','var')) 
        Domain = 1;
    end
    %-- default value for parameter thresh is 0
    if(~exist('thresh','var')) 
        thresh = 0;
    end  
    %-- default value for parameter color is 'r'
    if(~exist('color','var')) 
        color = 'r'; 
    end       
    %-- default behavior is to display intermediate outputs
    if(~exist('display','var'))
        display = true;
    end
    
    hfig = figure;
    if display
        axisRes = subplot(2,2,1);
        imshow(uint8(img));title('Contour evolution')
    else
       axisRes = subplot(1,1,1);
       imshow(uint8(img)); title('Contour evolution')
    end
    [row,col,channel] = size(img);
    
    nu = length*255*255; % coefficient of the length term
    
    c0=2;
    initialLSF = -init_mask.*c0 + (1 - init_mask).*c0;
    phi = initialLSF;
    OnesFilter = [];
    timestep = 0.2; % time step 
    mu = regularization; % coefficient of the level set (distance) regularization term P(\phi)
    epsilon = 5; % the paramater in the definition of smoothed Dirac function
    
    %% Select which domain to perform algorithm
    % Perform algorithm in Spatial Domain
    if Domain ==1 
           %set corresponding spatial domain filter
        if sigma2>=0.3
            filter2 = fspecial('gaussian',round(2*sigma2)*2+1,sigma2);% the Gaussian kernel 2  
        else
            filter2 = 1;
        end
        if sigma1>=0.3
            filter1 = fspecial('gaussian',round(2*sigma1)*2+1,sigma1);% the Gaussian kernel  1 
        else
            filter1 = 1;
        end
        OnesFilter{1} = conv2(ones(row,col),filter1,'same');
        OnesFilter{2} = conv2(ones(row,col),filter2,'same');
    % Perform algorithm in Frequency Domain   
    else
        D1= sigma1;D2=sigma2;
        [Fl1,Fl2]= freqspace([row,col],'meshgrid');   
        Rad = Fl1.^2+Fl2.^2;
        for i=1:row
           for j=1:col
                t1=Rad(i,j)/(2*D1*D1);t2=Rad(i,j)/(2*D2*D2);
                filter1(i,j)=exp(-t1);filter2(i,j)=exp(-t2);
           end
        end
%       [filter1  filter2] = CreatFilter(img,sigma1,sigma2,0.78);
      OnesFilter{1} =FreqFliter(ones(row,col),filter1);
      OnesFilter{2} =FreqFliter(ones(row,col),filter2);
     [flrow,flcol]    = size(filter1);filter = 1-filter2;
    end
   
    its = 0;      stop = 0;
    prev_mask = init_mask;  c = 0;
    Energy = [];energy = 10^9;
    while ((its < max_its) && ~stop)    
        %--
        if its==max_its-2
            its =its;
        end
        [phi energy,energy_img] = LSE_HFE(phi,img,nu,timestep,mu,Domain,...
            filter1,filter2,energy,epsilon,OnesFilter);
        Energy = [Energy,energy];
        new_mask = phi<=0;
        [c ndiff] = convergence(prev_mask,new_mask,thresh,c);
        if c <= 5
            its = its + 1;
            prev_mask = new_mask;
        else stop = 1;
        end      
        %-- intermediate output
        if ( mod(its,1)==0 )
            if display
                axisRes = subplot(2,2,1);
                showCurveAndPhi(phi,axisRes,color);
                subplot(2,2,2);
                surf(phi);colormap(summer)
                title('Level Set function')
                subplot(2,2,3),
                if its==1
                    imshow(energy_img,[]);
                    range_max = max(energy_img(:));
                    range_min = min(energy_img(:));
                else
                    imshow(energy_img,[range_min,range_max]);
                    range_max = max(energy_img(:));
                    range_min = min(energy_img(:));
                end
                subplot(2,2,4),
                imshow(double(phi<=0));title('Segmentation Results')
            else
                showCurveAndPhi(phi,axisRes,color);
            end
            drawnow;
            if (mod(its,60)==0)
                abc=its;
            end
        end
    end

    %-- final output
    if display
        showCurveAndPhi(phi,axisRes,color); 
        subplot(2,2,3),plot(Energy,'-kd','LineWidth',2,...
            'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',12);
        xlabel('Iteration','fontsize',12,'fontweight','b');
        ylabel('Energy','fontsize',12,'fontweight','b');title('Minimize FBE')
    else
        showCurveAndPhi(phi,axisRes,color);
    end
    its
    seg = phi<=0; %-- Get mask from levelset
    



%---------------------------------------------------------------------
%---------------------------------------------------------------------
%-- AUXILIARY FUNCTIONS ----------------------------------------------
%---------------------------------------------------------------------
%---------------------------------------------------------------------
  
%-- Displays the image with curve superimposed
function showCurveAndPhi(phi,axisRes,cl)

	axes(axisRes);
	delete(findobj(axisRes,'type','line'));
	hold on; [c,h] = contour(phi,[0 0],cl,'Linewidth',6);
%     [c1,h1] = contour(phi,[0 0],'w','Linewidth',3);
    hold off;
	delete(h);%delete(h1);
    test = isequal(size(c,2),0);
	while (test==false)
        s = c(2,1);%s1 = c1(2,1);
        if ( s == (size(c,2)-1) )
            t = c;%t1 = c1;
            hold on; plot(t(1,2:end)',t(2,2:end)',cl,'Linewidth',6);
            hold on;plot(t(1,2:end)',t(2,2:end)','w','Linewidth',3);
            test = true;
        else
            t = c(:,2:s+1);%t1 = c1(:,2:s1+1);
            hold on; plot(t(1,1:end)',t(2,1:end)',cl,'Linewidth',6);
            hold on;plot(t(1,1:end)',t(2,1:end)','w','Linewidth',3);
            c = c(:,s+2:end);
%             c1 = c1(:,s1+2:end);
        end
    end

% BFE: Author: Qing Guo
function [phi_new Energy_new,Energy_img] = LSE_HFE(phi_old,img,nu,timestep,mu,Domain,...
    filter1,filter2,Energy_old,epsilon,OnesFilter)
    [row,col,channel] = size(img);
    phi = phi_old;
    Energy = Energy_old+1;
    iter = 0;
    while (Energy-Energy_old>0)   
        phi = NeumannBoundCond(phi);
        K = curvature_central(phi);
        DrcU = (epsilon/pi)./(epsilon^2.+phi.^2);
        Hu = 0.5*(1+(2/pi)*atan(phi./epsilon));
        fb = zeros(row,col);
    for chi = 1: channel
        % Perform algorithm in spatial domain
        if Domain == 1
            [fb1_new,fb2_new] = bandpassitem(img(:,:,chi),phi,filter1,filter2,epsilon,OnesFilter);
            fb = -(fb1_new.^2.*Hu)+(fb2_new.^2.*(1-Hu))+fb;
        % Perform algorithm in Frequency domain
        else
            [fb1_new,fb2_new] = bandpassitemFre(img(:,:,chi),phi,filter1,filter2,epsilon,OnesFilter);
%             fb1_new_fft = fftshift(fft2(fb1_new)); fb2_new_fft = fftshift(fft2(fb2_new));
            fb = -(fb1_new.^2.*Hu)+(fb2_new.^2.*(1-Hu))+fb;
        end
    end
        item2 = (fb./channel).*DrcU;
        item3 = mu*(4*del2(phi)-K);
        item4 = nu.*DrcU.*K;
        phi = phi+timestep*(item2+item3+item4);
        Energy_img =(fb1_new.*Hu).^2+(fb2_new.*(1-Hu)).^2;
        Energy = sum(Energy_img(:));
        iter = iter +1;
        if iter > 1
            break;
        end
    end
    phi_new = phi;
    Energy_new = Energy;
    
%-- compute fh1 and fh2
function [fh1_new,fh2_new] = highpassitem(img,u,kernel,epsilon)
    Hu = 0.5*(1+(2/pi)*atan(u./epsilon));
    outside = img.*Hu;
    fh1_new = img-conv2(outside,kernel,'same')./conv2(Hu,kernel,'same');
    fh2_new = img-conv2(img.*(1-Hu),kernel,'same')./conv2(1-Hu,kernel,'same');
%-- compute fb1 and fb2
function [fb1_new,fb2_new] = bandpassitemFre(img,u,Hfilter1,Hfilter2,epsilon,OnesFilter)

    Hu = 0.5*(1+(2/pi)*atan(u./epsilon));  
    outside = img.*Hu; inside = img.*(1-Hu);
    lowoutside = FreqFliter(outside,Hfilter2);lowHuoutside = FreqFliter(Hu,Hfilter2);
    lowinside  = FreqFliter(inside,Hfilter2);lowHuinside = OnesFilter{2}-lowHuoutside;%FreqFliter(1-Hu,Hfilter2);
    highoutside = FreqFliter(outside,Hfilter1);highHuoutside = FreqFliter(Hu,Hfilter1);
    highinside  = FreqFliter(inside,Hfilter1);highHuinside = OnesFilter{1}-highHuoutside;%FreqFliter(1-Hu,Hfilter1);
    fl1_new = lowoutside./lowHuoutside;
    fl2_new = lowinside./lowHuinside;
    fh1_new = highoutside./highHuoutside;
    fh2_new = highinside./highHuinside;
    fb1_new = fh1_new-fl1_new;
    fb2_new = fh2_new-fl2_new;
    
function [fb1_new,fb2_new] = bandpassitem(img,u,kernel1,kernel2,epsilon,OnesFilter)
    Hu = 0.5*(1+(2/pi)*atan(u./epsilon));
    outside = img.*Hu;inside = img.*(1-Hu);
    HuKer1 = conv2(Hu,kernel1,'same');HuKer2 = conv2(Hu,kernel2,'same');
    fh1_new = conv2(outside,kernel1,'same')./HuKer1;
    fh2_new = conv2(inside,kernel1,'same')./(OnesFilter{1}-HuKer1);
    fl1_new = conv2(outside,kernel2,'same')./HuKer2;
    fl2_new = conv2(inside,kernel2,'same')./(OnesFilter{2}-HuKer2);
    fb1_new = fh1_new-fl1_new;
    fb2_new = fh2_new-fl2_new;
 
    
% Frequency domain filter
function [FreqRes] = FreqFliter(img,Hfilter)
  Fimg = fft2(img);
  Fimg = fftshift(Fimg);
  Fimg = Fimg.*Hfilter;
  Fimg = ifftshift(Fimg);
  FreqRes = real(ifft2(Fimg));
  
%-- Neumann boundary condition
function g = NeumannBoundCond(f)
    
    [nrow,ncol] = size(f);
    g = f;
    g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
    g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
    g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  
%-- compute curvature    
function k = curvature_central(u)                       

    [ux,uy] = gradient(u);                                  
    normDu = sqrt(ux.^2+uy.^2+1e-10);	% the norm of the gradient plus a small possitive number 
                                        % to avoid division by zero in the following computation.
    Nx = ux./normDu;                                       
    Ny = uy./normDu;
    nxx = gradient(Nx);                              
    [junk,nyy] = gradient(Ny);                              
    k = nxx+nyy;                        % compute divergence

% Convergence Test
function [c n_diff] = convergence(p_mask,n_mask,thresh,c)
    diff = p_mask - n_mask;
    n_diff = sum(abs(diff(:)));
    if n_diff < thresh
        c = c + 1;
    else c = 0;
    end

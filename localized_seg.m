%% Localizing Region Based Active Contour Segmentation:

function seg = localized_seg(I,init_mask,max_its,colorI,rad)
  
  %Weight of smoothing term
  alpha = .1; 
 
  display = true;
 
  I = double(I); 
 
  X = [];
  Y=[];
  %-- Default localization radius is 1/10 of average length
  [dimy dimx] = size(I);
  
    if(~exist('rad','var')) 
    rad = round((dimy+dimx)/(2*8)); 
    if(display>0) 
      disp(['localiztion radius is: ' num2str(rad) ' pixels']); 
    end
  end
   
  % Create a signed distance map (SDF) from mask
  phi = mask2phi(init_mask);  
  
  % take the narrow band as the initial contour
  idx = find(phi <= 0.5 & phi >= -0.5)';
    [dimy, dimx] = size(phi);  
  % save initial contour
  %########################################################change in
  %remaspgm
  initial_contour = zeros(dimy, dimx);
  initial_contour(idx) = 1;
  
   %  Calculate intial Energy
      int = find(phi<=0);
      int_avg = sum(I(int))/length(int);
      ext = find(phi>0);
      ext_avg = sum(I(ext))/length(ext);
      
      old_energy = (0.5 * (int_avg - ext_avg)^2);
      
  figure, imshow(logical(initial_contour));title('Initial Contour');
%   imwrite(initial_contour,'initial_contour.png');
%############################################################################3333  
  % do segmentation
  
  for its = 0:max_its   


    % get the curve's narrow band
    idx = find(phi <= 0.5 & phi >= -0.5)';
    [y x] = ind2sub(size(phi),idx);
      
    % get windows for localized statistics
    xneg = x-rad; xpos = x+rad;      
    yneg = y-rad; ypos = y+rad;
    xneg(xneg<1)=1; yneg(yneg<1)=1;  
    xpos(xpos>dimx)=dimx; ypos(ypos>dimy)=dimy;

    % re-initialize u,v,Ain,Aout
    u=zeros(size(idx)); v=zeros(size(idx)); 
    Ain=zeros(size(idx)); Aout=zeros(size(idx)); 
    
    % compute local stats
    for i = 1:numel(idx)  % for every point in the narrow band
      img = I(yneg(i):ypos(i),xneg(i):xpos(i)); %sub image
      P = phi(yneg(i):ypos(i),xneg(i):xpos(i)); %sub phi

      upts = find(P<=0);            %local interior
      Ain(i) = length(upts)+eps;
      u(i) = sum(img(upts))/Ain(i);
      
      vpts = find(P>0);             %local exterior
      Aout(i) = length(vpts)+eps;
      v(i) = sum(img(vpts))/Aout(i);
    end  
    
 
    % get image-based forces
    
    %F = (I(idx)-u).^2-(I(idx)-v).^2; %chan vase global
   %F = -(u-v).*(2.*I(idx)-u-v); %chan vase local
    
   F = -((u-v).*((I(idx)-u)./Ain+(I(idx)-v)./Aout)); % Yezzi local
    
    % get forces from curvature 
    curvature = get_curvature(phi,idx,x,y);  
    
    %gradient descent to minimize energy
    
    dphidt = F./max(abs(F)) + alpha*curvature;
    
    dt = .45/(max(dphidt)+eps);
    
    old_phi = phi;
    
    phi(idx) = phi(idx) + dt.*dphidt;
    
 %% Convergence using Stationarity of phi
 %{
    %difference in phi
    phi_diff = old_phi - phi;
    change = length(find(abs(phi_diff) >= 0.1)); 
    max_diff = max(max(phi_diff));
    
 %  set convergence threshold
    if(max_diff < 0.05)
        disp(['Converged in ', num2str(its),' iterations.']);
        break;
    end
 %}  
 %% Level set Reinitialization
    phi = sussman(phi, dt);
 
 %% Convergence using Energy minimization.
   
    %  Calculate New Energy
      int = find(phi<=0);
      int_avg = sum(I(int))/length(int);
      ext = find(phi>0);
      ext_avg = sum(I(ext))/length(ext);
      
      energy = -0.5 * (int_avg - ext_avg)^2;
      
      X = [X its];
      Y = [Y energy];
  
      change = old_energy - energy;     
      old_energy = energy;
      %{
      if change == 0 && its>1
        disp(['Converged in ', num2str(its),' iterations.']);
        break; 
      end
   %} 
  %% Display Intermediate Output
    
   
    % intermediate output
    if((display>0)&&(mod(its,50) == 0)) 
      showCurveAndPhi(I,phi,its);  
    end
   %}
    
  end % end of total iterations
 %% Display final output
  % final output
  
  if(display)
    showCurveAndPhi(I,phi,its);
  end
%}
  
 %% Plot energy Vs iterations
  % plot energy
   %plot(X,Y);
   
 %% Create final Segmentation Mask 
  
  % make mask from SDF
  seg = phi<=0; 

 

function showCurveAndPhi(img, phi, i)
  imshow(uint8(img)); hold on;
  contour(phi, [0 0], 'g','LineWidth',4);
  contour(phi, [0 0], 'k','LineWidth',2);
  hold off; title(['Segmentation Process: iteration ' num2str(i)]); drawnow;
  
% converts a mask to a SDF
function phi = mask2phi(init_a)
  phi=bwdist(init_a)-bwdist(1-init_a)+im2double(init_a)-.5;
  
% compute curvature along SDF
function curvature = get_curvature(phi,idx,x,y)
    
    [dimy, dimx] = size(phi);        

    ym1 = y-1; xm1 = x-1; yp1 = y+1; xp1 = x+1;

    ym1(ym1<1) = 1; xm1(xm1<1) = 1;              
    yp1(yp1>dimy)=dimy; xp1(xp1>dimx) = dimx;    

    idup = sub2ind(size(phi),yp1,x);    
    iddn = sub2ind(size(phi),ym1,x);
    idlt = sub2ind(size(phi),y,xm1);
    idrt = sub2ind(size(phi),y,xp1);
    idul = sub2ind(size(phi),yp1,xm1);
    idur = sub2ind(size(phi),yp1,xp1);
    iddl = sub2ind(size(phi),ym1,xm1);
    iddr = sub2ind(size(phi),ym1,xp1);
    
    % get central derivatives of SDF at x,y
    phi_x  = -phi(idlt)+phi(idrt);
    phi_y  = -phi(iddn)+phi(idup);
    phi_xx = phi(idlt)-2*phi(idx)+phi(idrt);
    phi_yy = phi(iddn)-2*phi(idx)+phi(idup);
    phi_xy = -0.25*phi(iddl)-0.25*phi(idur)...
             +0.25*phi(iddr)+0.25*phi(idul);
    phi_x2 = phi_x.^2;
    phi_y2 = phi_y.^2;
    
    % compute curvature 
    curvature = ((phi_x2.*phi_yy + phi_y2.*phi_xx - 2*phi_x.*phi_y.*phi_xy)./(phi_x2 + phi_y2 +eps).^(3/2)).*(phi_x2 + phi_y2).^(1/2);        
  

% level set re-initialization by the sussman method
function D = sussman(D, dt)
 
  a = D - shiftR(D); % backward
  b = shiftL(D) - D; % forward
  c = D - shiftD(D); % backward
  d = shiftU(D) - D; % forward
  
  a_p = a;  a_n = a; 
  b_p = b;  b_n = b;
  c_p = c;  c_n = c;
  d_p = d;  d_n = d;
  
  a_p(a < 0) = 0;
  a_n(a > 0) = 0;
  b_p(b < 0) = 0;
  b_n(b > 0) = 0;
  c_p(c < 0) = 0;
  c_n(c > 0) = 0;
  d_p(d < 0) = 0;
  d_n(d > 0) = 0;
  
  dD = zeros(size(D));
  D_neg_ind = find(D < 0);
  D_pos_ind = find(D > 0);
  dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
  dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;
  
  D = D - dt .* sussman_sign(D) .* dD;
  
% whole matrix derivatives
function shift = shiftD(M)
  shift = shiftR(M')';

function shift = shiftL(M)
  shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];

function shift = shiftR(M)
  shift = [ M(:,1) M(:,1:size(M,2)-1) ];

function shift = shiftU(M)
  shift = shiftL(M')';
  
function S = sussman_sign(D)
  S = D ./ sqrt(D.^2 + 1);    


	
	
I1=imread('/home/ramya/Desktop/Projects/Mitosis_cancer/weitest10.png');	
        %-- load the image
m = false(size(I1,1),size(I1,2));   %-- create initial mask
%m(37:213,89:227) = true;

%  m = false(size(I));          %-- create initial mask
  m(13:52,13:52) = true;
 
I1 = imresize(I1,.5);  %-- make image smaller 
m = imresize(m,.5);  %   for fast computation

subplot(2,2,1); imshow(I1); title('Input Image');
subplot(2,2,2); imshow(m); title('Initialization');
subplot(2,2,3); title('Segmentation');

seg = localized_seg(I1, m, 400);  %-- run segmentation

subplot(2,2,4); imshow(seg); title('Final Segmentation');
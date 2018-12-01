function threshold
out='/home/ramya/Desktop/Projects/Mitosis_cancer';
I=imread('/home/ramya/Desktop/Projects/Mitosis_cancer/output1/A03_00Bb8true.png');
[threshold I1]=maxentropie(I)
[m n]=size(I1);

      
            
baseFileName = sprintf('%s.png','test10' );
    fullFileName = fullfile(out, baseFileName); 
    imwrite(I1, fullFileName);
imshowpair(I,I1,'montage');
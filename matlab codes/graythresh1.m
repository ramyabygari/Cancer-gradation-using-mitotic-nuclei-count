I=imread('/home/ramya/Desktop/Projects/Mitosis_cancer/redoutput1/A03_00Bb8true.png');
level=graythresh(I)
BW=imbinarize(I,level)
imshowpair(I,BW,'montage')
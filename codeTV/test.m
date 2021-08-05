g1=load('phanton_meas.mat');
g2=g1.meas;
%g3= uint8(255 * g2);
%map = hsv(256); 
%g4 = ind2rgb(g3, map); 
dTV=gradTVcc(g2);
dl2=gradlp(g2,2);
%dlhalf = gradlp(g2,0.5);

subplot(3,1,1)
imshow(g2)
subplot(3,1,2)
imshow(dTV)
subplot(3,1,3)
imshow(100*dl2)
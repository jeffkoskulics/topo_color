x = z(:,:,1);
%{
k =2*pi/0.0505;
%}
A = 0.005;
phi = -0.005;
f = @(x) A*sin(k*(x-phi));
sine_profile = f(x);
sine_profile(1:99,:) = 0;
sine_profile(559:1025,:) = 0;
sine_profile(:,1:105)= 0;
sine_profile(:,968:1281) = 0;

R = sine_profile - z(:,:,3);
figure; imagesc(R);
hold off
figure; plot(z(500,:,1),sine_profile(500,:),' - ','Color','Blue')
hold on
plot(z(200,:,1),z(150,:,3),'Color','Black')
%plot(z(250,:,1),z(250,:,3),'Color','Green')
%plot(z(350,:,1),z(350,:,3),'Color','Cyan')
%plot(z(450,:,1),z(450,:,3),'Color','Magenta')
%plot(z(500,:,1),z(500,:,3),'Color','Cyan')
set(gca,'DataAspectRatio',[1 1 1])
fid = fopen('s.out');
stuff = fscanf(fid,'%g',[5,Inf]);
plot(stuff(1,:),stuff(2,:),'1');
xlabel('t')
ylabel('a')
hold on
plot(stuff(1,:),stuff(3,:),'2');
plot(stuff(1,:),stuff(4,:),'3');
plot(stuff(1,:),stuff(5,:),'4');
legend('a11','a22','a33','a12')
hold off
print('s-out.eps','-deps');

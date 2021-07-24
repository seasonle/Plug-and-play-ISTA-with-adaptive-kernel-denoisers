


Spsnrs1=Spsnrs(:,1:300);
Dpsnrs1=Dpsnrs(:,1:300);


% Dpsnrs1(:,1:10)=Spsnrs(:,1:10);
% 
% 
% Dpsnrs1(:,11:100)=Dpsnrs(:,11:100);
figure;
  plot(Spsnrs1(1,:),'--r','LineWidth',1)%%0.3

  hold on;
 plot(Spsnrs1(2,:),'--b','LineWidth',1)%%0.9

  hold on;
 plot(Spsnrs1(3,:),'--g','LineWidth',1)%%1.2
 hold on;
 
 plot(Spsnrs1(4,:),'--k','LineWidth',1)%%1.5
 hold on;
 plot(Spsnrs1(5,:),'--m','LineWidth',1)%%1.5
 

 hold on;
 plot(Dpsnrs1(1,:),'r','LineWidth',1)%%0.3
  hold on;
 plot(Dpsnrs1(2,:),'b','LineWidth',1)%%0.9
  hold on;
 plot(Dpsnrs1(3,:),'g','LineWidth',1)%%1.2
  hold on;
 plot(Dpsnrs1(4,:),'k','LineWidth',1)%%1.5
 hold on;
 plot(Dpsnrs1(5,:),'m','LineWidth',1)%%1.5
 
grid on; axis tight; xlabel('Iteration'); ylabel('PSNR');
% title('PSNR');
 legend('r=0.6(F)','r=0.9(F)','r=1.2(F)','r=1.6(F)','r=1.8(F)','r=0.6(A)','r=0.9(A)','r=1.2(A)','r=1.6(A)','r=1.8(A)')
% legend('r=0.3(F)','r=0.5(F)','r=0.7(F)','r=0.9(F)','r=1.2(F)','r=0.3(A)','r=0.5(A)','r=0.7(A)','r=0.9(A)','r=1.2(A)')


% legend('r=0.3(F)','r=0.5(F)','r=0.7(F)','r=0.9(F)','r=0.3(A)','r=0.5(A)','r=0.7(A)','r=0.9(A)')

% legend('r=0.3(F)','r=0.5(F)','r=0.9(F)','r=1.9(F)','r=0.3(A)','r=0.5(A)','r=0.9(A)','r=1.9(A)')
% legend('r=0.3(F)','r=0.5(F)','r=0.9(F)','r=1.5(F)','r=0.3(A)','r=0.5(A)','r=0.9(A)','r=1.5(A)')

% legend('r=0.3(F)','r=0.9(F)','r=1.5(F)','r=0.3(A)','r=0.9(A)','r=1.5(A)')

% legend('r=0.3(F)','r=0.9(F)','r=0.3(A)','r=0.9(A)')

% axis ( [0 100 24.3 28] );

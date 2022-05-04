clc
clear all
load('data_gen.mat')
name='boat'
x0 = load_image(name,M);
x0 = rescale( sum(x0,3) );
x=zeros(3*M*N,1);
x(1:M*N)=reshape(x0,M*N,1);
b=A*x;
clf;
subplot(1,2,2)
imageplot(reshape(b(1:M*N),M,N));
subplot(1,2,1)
imageplot(x0)
norm(x0-reshape(b(1:M*N),M,N),'fro')
save(name,"b","x0","x")
% y=x0;
% a=0;
% bb=1;
% nn=200;
% D=a*eye(N,N);
% for i=1:(N-1)
%     D(i,i+1)=bb;
%     D(i+1,i)=bb;
% end
% D1=(1/(a+3*bb))*eye(M,M);
% D2=(1/(a+4*bb))*eye(M,M);
% D1(1,1)=1/(a+2*bb);
% D2(1,1)=1/(a+3*bb);
% D1(M,M)=1/(a+2*bb);
% D2(M,M)=1/(a+3*bb);
% K1=D1*D;
% K2=bb*D1;
% K3=D2*D;
% K4=bb*D2;
% for kkk=1:10
%     y_aux=y;
%  for i=1:(N)
%     if i==1
%       y(:,i)=K1*y_aux(:,1)+K2*y_aux(:,2);  
%     elseif i==M
%         y(:,i)=K1*y_aux(:,M)+K2*y_aux(:,M-1);
%     else
%         y(:,i)=K3*y_aux(:,i)+K4*y_aux(:,i-1)+K4*y_aux(:,i+1);
%     end
%  end
% end


clear all; clc;

nn=200;
M=2*nn; N=5*nn;
K=rand(M,N);
x_gen=linspace(-1,1,1000);
x_gen=atan(x_gen)';
clf
plot(x_gen)
D=zeros(N,N)
D(1,1)=1;
D(N,1)=1;
for i=1:(N-1)
    D(i,i+1)=1;
    D(i+1,i)=-1;
end
A=[K zeros(M,N); D -eye(N,N)]
g=K*x_gen;
b=[g;zeros(N,1)]
% coherence=0;
% for i = 1 : M-1
%     for j = i+1 : M
%         ccc = abs(A(:,i)'*A(:,j));
%         if ccc > coherence
%             coherence = ccc;
%         end
%     end
% end
% 
% coherence
% 
% % 2s-1 < 1/coher
% % s<(1+1/coher)/2
% 
% % sgen=(coherence+1)/(2*coherence)
% sgen=N/20
% x_gen=sprand(N,1,sgen/N);
% s=nnz(x_gen)
% 2*s-1 < 1/coherence
% b=A*x_gen;
% 
% 
% 
% v_cost=[zeros(1,N),ones(1,N)];
% C_P=[-eye(N), -eye(N); eye(N), -eye(N)];
% b_P=zeros(2*N,1);
% Aeq=[A, zeros(M,N)];
% [u, val_P] = linprog(v_cost, C_P, b_P, Aeq, b);
% x_exact=u(1:N,1); 
% 
% norm(x_exact-x_gen)
% norm(pinv(A)*b-x_gen)
% 
% C_D = [A'; -A'];
% b_D=ones(2*N,1);
% [mu_exact, val_D] = linprog(b, C_D, b_D);
% 
% 
% norm(A*x_exact-b)
% 
% w=-A'*mu_exact;
% toll=10^-7;
% count=0;
% for i=1:N
%     if (x_exact(i)<-toll && norm(w(i)+1)>toll)
%         count=count+1;
%     elseif (x_exact(i)>toll && norm(w(i)-1)>toll)
%         count=count+1;
%     elseif (norm(x_exact(i))<=toll && ( w(i)>=1+toll || w(i)<=-1-toll))
%         count=count+1;
%     end
% end
% count
% 
% 
% save('data_gauss')
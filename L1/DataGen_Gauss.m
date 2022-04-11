clear all; clc;
rng();%Same seed
nn=600;
tic
M=2*nn; N=5*nn;%Matrix size
M=2260; N=3000;%Matrix size
A=zeros(M,N);
% Matrix such that every entry is an independent sample from normal
% Then is normalized column by column
for i=1:M
    for j=1:N
        A(i,j)=normrnd(0,1);
    end
end
 for j=1:N
     A(:,j)=A(:,j)/norm(A(:,j));
end
sgen=N/10%Sparsity
x_gen=sprand(N,1,sgen/N);
b=A*x_gen;%data



v_cost=[zeros(1,N),ones(1,N)];
C_P=[-eye(N), -eye(N); eye(N), -eye(N)];
b_P=zeros(2*N,1);
Aeq=[A, zeros(M,N)];
[u, val_P] = linprog(v_cost, C_P, b_P, Aeq, b);%Linear problem
x_exact=u(1:N,1); 

norm(x_exact-x_gen)
norm(pinv(A)*b-x_gen)

C_D = [A'; -A'];
t=1
b_D=ones(2*N,1);
[mu_exact, val_D] = linprog(b, C_D, b_D);%Linear problem parameters
t-2

norm(A*x_exact-b)

w=-A'*mu_exact;
toll=10^-14;
count=0;
for i=1:N%Number of entry satifisfying the first order condition
    if (x_exact(i)<-toll && norm(w(i)+1)>toll)
        count=count+1;
    elseif (x_exact(i)>toll && norm(w(i)-1)>toll)
        count=count+1;
    elseif (norm(x_exact(i))<=toll && ( w(i)>=1+toll || w(i)<=-1-toll))
        count=count+1;
    end
end
count
toc
save('data_gauss')% Saving the data
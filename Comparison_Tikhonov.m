clc 
clear all
load('data_gauss');
rng(0);
MaxIt = 300;% Maximun number of Iteration
normA=norm(A);%Norm of A
delta=0.4*norm(b);%Level of Noise
db=-0.5+rand(M,1); db=db/norm(db);
tol=10^-3;%tolerance
b_delta=b+delta*db;%Noisy data
x_exact=x_gen;
penalty=0;%Penalty term
xk=zeros(N,1); x_old=xk;%Initialization
lambda=0.999/normA^2;%Step sizes
deep=6;
Length=5;%Grid number
ES=zeros(N,1);
Penalties=zeros(1,Length*deep+1);%Sequences of different penalties
Errors=zeros(1,Length*deep+1);%Sequences of error obtained with diferent penalties;
solutions=zeros(N,Length*deep+1);%Sequences of solutions(vector) obtained with diferent penalties;
tic
for d=1:deep 
    for l=1:Length
        penalty=(1-(l-1)/Length)*10^(1-d);
        xk=ES;
        x_old=xk;
        % Tikhonov FB
        %Dist: Norm of the iteration with respect the solution
        %penalty: Penalty term
        %EST: Execution time
        %ESI: Number of iteration
        %error: Mean square error
        [Dist,penalty,EST,ESI,error]=Tikhonov(xk,lambda,MaxIt,A,b_delta,b,x_exact,normA,tol,penalty,M,N,delta);
        Errors(1,((d-1)*Length+l))=min(Dist); %Errro
        Penalties(1,((d-1)*Length+l))=penalty;
        solutions(:,((d-1)*Length+l))=ES;
        (d-1)*Length+l
    end
end
[Dist,penalty,ES,EST,ESI,error] =Tikhonov(xk,lambda,MaxIt,A,b_delta,b,x_exact,normA,tol,0,M,N,delta);
Errors(1,Length*deep+1)=min(Dist);%Sequence of errors
Penalties(1,Length*deep+1)=0;%Sequences of penalties
solutions(:,Length*deep+1)=ES;%Sequence of solutions
toc%time
clf%Figure
figure; semilogx(Penalties, Errors,'b')
legend('Tikhonov','Interpreter','latex')
title(strcat('M =  ', num2str(M),', N =  ', num2str(N), 'Interpreter','latex'))
xlabel('Penalty (\lambda)') 
ylabel('$\|x_{\lambda}-x_*\|$','Interpreter','latex') 
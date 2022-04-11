
function [BD,Feas,Dist,ValP,EST,ESI,Time_Tol,Iter_Tol,Time,error]= DouglRach_primal(xk,lambda,MaxIt,A,b_delta,b,x_exact,mu,normA,tol)
ESI=0;% Iteration number
EST=0;% Execution time
Time_Tol=0; %Time until reach the tolerance
Iter_Tol=0; %Iteration until reach the tolerance
Time=0;%Total time
error=norm(xk-x_exact);
BD=zeros(MaxIt,1);%Bregman distance % BD_ERG=zeros(MaxIt,1);
Feas=zeros(MaxIt,1);%Feasibility % Feas_ERG=zeros(MaxIt,1);
Dist=zeros(MaxIt,1);%Norm error % Dist_ERG=zeros(MaxIt,1);
ValP=zeros(MaxIt,1);%Function value

N=length(A(1,:)); yk=xk;%Initialization
% PseudoAt=pinv(A');
PseudoA=pinv(A); %Pseudo inverse
tic%
% Pseudo_b=PseudoA*b_delta;
for k=1:MaxIt
        for i=1:N
            if yk(i,1)<-lambda
                xk(i,1)=yk(i)+lambda;
            elseif yk(i,1)>lambda
                xk(i,1)=yk(i)-lambda;
            else
                xk(i,1)=0;
            end
        end
    
        x_aux=2*xk-yk;

        p=x_aux-PseudoA*(A*x_aux-b_delta);
        yk=yk+p-xk;
        if norm(yk-x_exact)<error
            EST=toc;
            ESI=k;
            error=norm(yk-x_exact);
        end
        if norm(p-x_exact)<error
            EST=toc;
            ESI=k;
            error=norm(p-x_exact);
        end
        if norm(xk-x_exact)<error
            EST=toc;
            ESI=k;
            error=norm(xk-x_exact);
        end
        if norm(xk-x_exact)>tol
            Time_Tol=toc;
            Iter_Tol=k;
        end
        ValP(k,1)=abs(sum(abs(xk))-sum(abs(x_exact)));%Function value
        BD(k,1)=abs(sum(abs(xk))-sum(abs(x_exact))+ mu'*(A*xk-b_delta)); %Bregman distance
        Feas(k,1)=norm(A*xk-b); % Feasibility
        Dist(k,1)=norm(xk-x_exact); % Reconstruction error
        
end
Time=toc;
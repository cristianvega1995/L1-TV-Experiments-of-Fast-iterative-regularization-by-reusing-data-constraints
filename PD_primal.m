
function [BD,Feas,Dist,ValP,EST,ESI,Time_Tol,Iter_Tol,Time,error] = PD_primal(xk,xold,uk,lambda,MaxIt,A,b_delta,b,x_exact,mu,normA,tol)
ESI=0;% Iteration number
EST=0;% Execution time
Time_Tol=0; %Time until reach the tolerance
Iter_Tol=0; %Iteration until reach the tolerance
Time=0;%Total time
error=norm(xk-x_exact);%Error
BD=zeros(MaxIt,1);%Bregman distance % BD_ERG=zeros(MaxIt,1);
Feas=zeros(MaxIt,1);%Feasibility % Feas_ERG=zeros(MaxIt,1);
Dist=zeros(MaxIt,1);%Norm error % Dist_ERG=zeros(MaxIt,1);
ValP=zeros(MaxIt,1);%Function value

N=length(A(1,:));
tic
for k=1:MaxIt
        uk=uk+lambda*(A*(2*xk-xold)-b_delta);
        
        xold=xk;
        x_aux=xk-lambda*A'*uk;
        for i=1:N
            if x_aux(i,1)<-lambda
                xk(i,1)=x_aux(i)+lambda;
            elseif x_aux(i,1)>lambda
                xk(i,1)=x_aux(i)-lambda;
            else
                xk(i,1)=0;
            end
        end
        if norm(xk-x_exact)>tol
            Time_Tol=toc;
            Iter_Tol=k;
        end
        if norm(xk-x_exact)<error
            EST=toc;
            ESI=k;
            error=norm(xk-x_exact);
        end
        ValP(k,1)=abs(sum(abs(xk))-sum(abs(x_exact)));%Function value
        BD(k,1)=abs(sum(abs(xk))-sum(abs(x_exact))+ mu'*(A*xk-b_delta)); %Bregman distance
        Feas(k,1)=norm(A*xk-b); % Feasibility
        Dist(k,1)=norm(xk-x_exact); % Reconstruction error  
end
Time=toc;
 
function [Dist,penalty,EST,ESI,error] =Tikhonov(xk,lambda,MaxIt,A,b_delta,b,x_exact,normA,tol,penalty,M,N,delta)
ESI=0;% Iteration number
EST=0;% Execution time
error=norm(xk-x_exact);
Feas=zeros(MaxIt,1);%Feasibility
Dist=zeros(MaxIt,1);%Norm error
ValP=zeros(MaxIt,1);%Function value
penalty=penalty*max(abs(A'*b_delta));%Penalty

lambda=0.99/normA^2;%Step size
N=length(A(1,:));
tic
k=0;
error1=error;
xold=xk;
while error1>tol & k<MaxIt
 k=k+1;       
        x_aux=xk-lambda*A'*(A*xk-b_delta);
        for i=1:N
            if x_aux(i,1)<-lambda*penalty
                xk(i,1)=x_aux(i)+lambda*penalty;
            elseif x_aux(i,1)>lambda*penalty
                xk(i,1)=x_aux(i)-lambda*penalty;
            else
                xk(i,1)=0;
            end
        end
        if norm(xk-x_exact)<error
            EST=toc;
            ESI=k;
%             ES=xk;
            error=norm(xk-x_exact);
        end
        error1=norm(xk-xold);
        xold=xk;
        Dist(k,1)=norm(xk-x_exact); % Dist_erg(k,1)=norm(Xk-x);
        
end
end
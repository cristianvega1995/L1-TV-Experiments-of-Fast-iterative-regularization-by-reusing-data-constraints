function [BD,Feas,Dist,ValP,EST,ESI,Time_Tol,Iter_Tol,Time,error]= PROJ_primalseries(xk,pk,p_old,uk,lambda,prob,MaxIt,A,b_delta,b,x_exact,mu,normA,tol)
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

N=length(A(1,:)); M=length(A(:,1));%Matrix size

p_bar=pk+xk-p_old; % Initialization
tic;
for k=1:MaxIt
        
    
        uk=uk+lambda*(A*p_bar-b_delta);
        
        x_aux=pk-lambda*A'*uk;
        for i=1:N
            if x_aux(i,1)<-lambda
                xk(i,1)=x_aux(i)+lambda;
            elseif x_aux(i,1)>lambda
                xk(i,1)=x_aux(i)-lambda;
            else
                xk(i,1)=0;
            end
        end
        
        
        p_old=pk;
        pk=xk;
        mm=randperm(M);%Permutation
        for j=1:M%Composition of projections
           
            if rand < prob
               
                indx=mm(j);
                a=A(indx,:);
                pk = pk-  ( (a*pk-b_delta(indx))/(norm(a)^2) ) * a';%Projections
                
            end
            
        end
        
        p_bar=pk+xk-p_old;
        if norm(xk-x_exact)>tol
            Time_Tol=toc;

        end
        if norm(xk-x_exact)>tol%Tolerance
            Time_Tol=toc;
            Iter_Tol=k;
        end
        if norm(xk-x_exact)<error%Stopping time
            EST=toc;
            ESI=k;
            error=norm(xk-x_exact);
        end
        ValP(k,1)=abs(sum(abs(xk))-sum(abs(x_exact)));%Function value
        BD(k,1)=abs(sum(abs(xk))-sum(abs(x_exact))+ mu'*(A*xk-b_delta)); %Bregman distance
        Feas(k,1)=norm(A*xk-b); % Feasibility
        Dist(k,1)=norm(xk-x_exact); % Reconstruction error
        
end
Time=toc;%Execution time
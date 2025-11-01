      %parameters
      
Nmax=10; %Number of lattice sites
NTest=1000;%Number of simulations
JJ=40; %An internal parameter that allows us to study timeascale values
ts=0; %Time initial condition
c = -1/2; %coupling parameter
h = 50; %Inertia parameter 
f=1; % Additional parameter to study the Landau transition
s = 1; %wienner process size

n=Nmax;
Defno4=zeros(NTest,JJ); %Size of exported table


%Parallelisation of the simulations for different timescales

parfor jl=1:JJ
    
    tq = 10^(jl/10);
    dt = 0.01;
    Ns= floor(tq/dt);
    for z=1:NTest
    
%Random numbers
        Wn=zeros(n,Ns+1);
        for k=1:n
            for i=1:Ns+1
                Wn(k,i)=randn();
            end
        end
%Creating the initial condition

        Xt=zeros(n,Ns);
        Pt=zeros(n,Ns);

%Solving the equation
        t=ts;

        for i=2:Ns
            Xs=zeros(n,1);
            Ps=zeros(n,1);
            for k=1:n
                if k==1
                    Ps(1)=Pt(1,i-1)-(h*Pt(1,i-1)+c*Xt(2,i-1)+f*(Xt(1,i-1)*(-2*(t./tq)+1)+Xt(1,i-1)^3) )*dt + s*Wn(1,i-1)*sqrt(dt);
                    Xs(1)=Xt(1,i-1)+Pt(1,i-1)*dt;
                else
                    if k==n
                        Ps(n)=Pt(n,i-1)-(h*Pt(n,i-1)+c*Xt(n-1,i-1)+f*(Xt(n,i-1)*(-2*(t./tq)+1)+Xt(n,i-1)^3) )*dt + s*Wn(n,i-1)*sqrt(dt);
                        Xs(n)=Xt(n,i-1)+Pt(n,i-1)*dt;
                    else
                        Ps(k)=Pt(k,i-1)-(h*Pt(k,i-1)+c*(Xt(k-1,i-1)+Xt(k+1,i-1))+f*(Xt(k,i-1)*(-2*(t./tq)+1)+Xt(k,i-1)^3) )*dt + s*Wn(k,i-1)*sqrt(dt);
                        Xs(k)=Xt(k,i-1)+Pt(k,i-1)*dt;
                    end
            
                end
        
            end
            Xt(:,i)=Xs;
            Pt(:,i)=Ps;
            t=t+dt;
        end

        defno=0;
        for j=1:n-1
            if Xt(j,Ns)*Xt(j+1,Ns)<0
                defno=defno+1;
            end
        end
        Defno4(z,jl)=defno;
    end
end


writematrix(Defno4,"Defno_comp10.txt");
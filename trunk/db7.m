%IMM based multiuser detector for Synchronous channel with two transmit antennae per user and two recieve antennae system.
%with Alamouti's space-time coding, no multipath, time-varying channel gains fixed for all simulations, slow, flat Rayleigh 
%fading channel (coherence BW > channel BW). Training based.
%Date : 5 March

clear all;
%format long e;
K = 3;                      %no of users
P = 2;                      %receiver antennae
n = 2;
T = 50000;                    %length of user bit stream
A = 5*eye([K*4 K*4]);       %amplitude of the kth users signal 
Q   = 4e-3 * eye(4);
db =20;   %-2
L=1;

gold_seq(1,:,1)= [-1 -1 -1 -1 1 1 -1 1 1 -1 -1 1 -1 1 -1 -1 1 1  1 1 1 -1 1 1 -1 -1 1 -1 -1 1 1];            % 1 transmitter for user 1
gold_seq(1,:,2)= [1 -1 1 -1 -1 -1 -1 1 1 -1 1 -1 1 1 -1  1 -1 1 -1 -1 -1 1 -1 1 1 -1 1 -1 1 1 1];            % 2 transmitter for user 1
gold_seq(2,:,1)= [-1 -1 1 -1 1 -1 -1 -1 1 1 1 1 -1 1 -1 1 -1 -1 1 1 -1 1 1 -1 1 1 -1 1 -1 -1 1];
gold_seq(2,:,2)= [-1 1 1 -1 1 1 -1 -1 -1 1 -1 1 1 -1 -1 1 -1 -1 -1 -1 1 1 1 1 -1 1 1 -1 1 1 -1];
% gold_seq(2,:,1)=[-1     1    -1     1     1    -1     1     1     1    -1    -1    -1     1  -1    -1    -1    -1    -1     1    -1    -1     1    -1    -1    -1     1   1    -1    -1    -1     1];
% gold_seq(2,:,2)=[-1     1     1     1    -1    -1    -1    -1     1    -1    -1    -1    -1  -1     1     1     1    -1    -1    -1    -1     1    -1    -1    -1    -1   -1    -1    -1    -1    -1];
gold_seq(3,:,1)= [-1 1 1 -1 -1 -1 1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 1 -1 1 1 -1 1 -1 -1 1 1 1 -1 1];
gold_seq(3,:,2)= [-1 1 -1 -1 -1 1 1 1 -1 1 -1 1 -1 1 1 1 -1 1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 1 1];

Tc=31;
fc=1/Tc;
nv=0;
gold_seq=gold_seq/sqrt(Tc);

[Bf,Af]=butter(3,0.03);                         %3rd order butterworth filter with fading rate 2*0.07
i=sqrt(-1);
G=randn(2*T+400,K*4*2*L);                                  % +1 for Rician
G=filter(Bf,Af,G);
g(1,:)=(G(:,1)-G(:,2)*i)';
g(2,:)=(G(:,3)-G(:,4)*i)';
g(3,:)=(G(:,5)-G(:,6)*i)';
g(4,:)=(G(:,7)-G(:,8)*i)';
g(5,:)=(G(:,9)-G(:,10)*i)';
g(6,:)=(G(:,11)-G(:,12)*i)';
g(7,:)=(G(:,13)-G(:,14)*i)';
g(8,:)=(G(:,15)-G(:,16)*i)';
g(9,:)=(G(:,17)-G(:,18)*i)';
g(10,:)=(G(:,19)-G(:,20)*i)';
g(11,:)=(G(:,21)-G(:,22)*i)';
g(12,:)=(G(:,23)-G(:,24)*i)';

%input bit stream for each user
b(:,:)=randint1(2*K,T);         
b(:,:)=b(:,:)-(b(:,:)==0);

%multipath delay tau(k,1,q,p) - the rx antenna shld be spaced sufficiently apart so that multipath components have sig. diff. prop. delays. 
tau(1,1,:,:) = [0 0; 0 0];              % row - transmitter col-  receiver
tau(2,1,:,:) = [0 0; 0 0];
tau(3,1,:,:) = [0 0; 0 0];


%time-varying channel gain for kth user lth path, mo multipath


%sig seq matrix S1 for time i-1
for p=1:P,
    for k1=1:K,
        for q1=1:2,
            for p1=1:P,
                if (p == p1)
                    d=tau(k1,1,q1,p1);
                    %S1(p, (k1-1)*4*Tc + (q1-1)*2*Tc + (p1-1)*Tc +1: (k1-1)*4*Tc + (q1-1)*2*Tc + p1*Tc) = [gold_seq(k1,(d+1):Tc,q1) zeros(1,d)];
                    S1(p, (k1-1)*4*Tc + (q1-1)*2*Tc + (p1-1)*Tc +1: (k1-1)*4*Tc + (q1-1)*2*Tc + p1*Tc) = [ gold_seq(k1,(Tc-d)+1:Tc,q1) gold_seq(k1,1:(Tc-d),q1)]; % nel correlation between transmitter antennae
                else
                    S1(p, (k1-1)*4*Tc + (q1-1)*2*Tc + (p1-1)*Tc +1: (k1-1)*4*Tc + (q1-1)*2*Tc + p1*Tc) = zeros(1,Tc);       %As nel correlation between rx antennae is assumed zero
                end;
            end;
        end;
    end;
end;    
        
S2=S1';
for r=1:K*4,
    for c1=1:K*4,
        W(r,c1)=0;
        for c=1:2,
            W(r,c1)=W(r,c1) + S1(c,(c1-1)*Tc + 1 : c1 * Tc) * S2((r-1)*Tc + 1 : r * Tc,c);
        end;
    end;
end;

%Matrix factorisation UL
UP=ULfact(2*K*2,W);
LT = UP';     

%sig seq matrix S3 for time i
for p=1:P,
    for k1=1:K,
        trn=2;
        for q1=1:2,
            for p1=1:P,
                if (p == p1)
                    d=tau(k1,1,q1,p1);
                    %S3(p, (k1-1)*4*Tc + (q1-1)*2*Tc + (p1-1)*Tc +1: (k1-1)*4*Tc + (q1-1)*2*Tc + p1*Tc) = [gold_seq(k1,(d+1):Tc,trn) zeros(1,d)];
                    S3(p, (k1-1)*4*Tc + (q1-1)*2*Tc + (p1-1)*Tc +1: (k1-1)*4*Tc + (q1-1)*2*Tc + p1*Tc) = [gold_seq(k1,(Tc-d)+1:Tc,trn) gold_seq(k1,1:(Tc-d),trn) ];

                else
                    S3(p, (k1-1)*4*Tc + (q1-1)*2*Tc + (p1-1)*Tc +1: (k1-1)*4*Tc + (q1-1)*2*Tc + p1*Tc) = zeros(1,Tc);
                end;
            end;
            trn=trn-1;
        end;
    end;
end;    
S4=S3';
for r=1:K*4,
    for c1=1:K*4,
        W1(r,c1)=0;
        for c=1:2,
            W1(r,c1)=W1(r,c1) + S3(c,(c1-1)*Tc + 1 : c1 * Tc) * S4((r-1)*Tc + 1 : r * Tc,c);
        end;
    end;
end;
%Matrix factorisation UL
UP1=ULfact(2*K*2,W1);
LT1 = UP1';
%sig1=A(1,1)^2*( W(1,1)^2*g(1,:)*ctranspose(g(1,:)) + W(3,3)^2*g(3,:)*ctranspose(g(3,:))+ W(2,2)^2*g(2,:)*ctranspose(g(2,:)) + W(4,4)^2*g(4,:)*ctranspose(g(4,:)) );
% for db=-2:2:7,
    nv=nv+1;
    V=10^(-db/10);                              %Additive noise variance
    Rchmse = zeros(4,2*T);
    Ichmse = zeros(4,2*T);
    
    for ppp=1:20,  %change parameters below                              %Simulations are repeated 20 times 
        
        
        v=sqrt(V)*randn(T,2,1);		            %Observation noise at time t-1
        v=v(:,1,:)+v(:,2,:)*sqrt(-1);           %Complex noise variance is 2V
        v=squeeze(v); 
        v1(:,1)=v(:,1);
        
        v=sqrt(V)*randn(T,2,1);		            %Observation noise at time t-1
        v=v(:,1,:)+v(:,2,:)*sqrt(-1);
        v=squeeze(v); 
        v1(:,2)=v(:,1);
        
        v=sqrt(V)*randn(T,2,1);		            %Observation noise at time t-1
        v=v(:,1,:)+v(:,2,:)*sqrt(-1);
        v=squeeze(v); 
        v2(:,1)=v(:,1);
        
        v=sqrt(V)*randn(T,2,1);		            %Observation noise at time t-1
        v=v(:,1,:)+v(:,2,:)*sqrt(-1);
        v=squeeze(v); 
        v2(:,2)=v(:,1);
        
        %Observation noise at time t
    
        for t=1:T,
            %for time t-1
            Bt = diag([kron(b(:,t),ones(2,1))]);
            tmp = A * Bt * g(:,(t-1)*2+1+400);
            nt  = [v1(t,1) * ones(Tc,1); v1(t,2) * ones(Tc,1)];
            n   = zeros(K*4,1);
            
            for q=1:K*4,
                n(q) = [S2((q-1)*Tc +1 : q*Tc,1)' S2((q-1)*Tc +1 : q*Tc,2)' ] * nt;
            end;
            
            z(:,t)=W * tmp + n;
            %noise whitening
            y1(:,t) = inv(UP)* z(:,t);
            
            %for time t
%             Bt=blkdiag(-b(2,t)*eye(2) , b(1,t)*eye(2) , -b(4,t)*eye(2) , b(3,t)*eye(2) , -b(6,t)*eye(2), b(5,t)*eye(2));
            Bt=blkdiag(-b(2,t)*eye(2) , b(1,t)*eye(2) , -b(4,t)*eye(2) , b(3,t)*eye(2), -b(6,t)*eye(2) , b(5,t)*eye(2));
            tmp = A * Bt * g(:,t*2+400);
            nt  = [v2(t,1) * ones(Tc,1); v2(t,2) * ones(Tc,1)];
            n   = zeros(K*4,1);
            
            for q=1:K*4,
                n(q)   =  [S4((q-1)*Tc +1 : q*Tc,1)' S4((q-1)*Tc +1 : q*Tc,2)' ]* nt;
            end;

            z1(:,t)=W1 * tmp + n;

            %noise whitening
            y2(:,t) = inv(UP1)* z1(:,t);
            
        end;
        snr = (A(1,1)^2*( g(1,401:T+400)*ctranspose(g(1,401:T+400)) + g(2,401:T+400)*ctranspose(g(2,401:T+400)) + g(3,401:T+400)*ctranspose(g(3,401:T+400)) + g(4,401:T+400)*ctranspose(g(4,401:T+400)) ) )/(4*T*2*V);
        SNR(ppp,nv) = 10*log10(real(snr));
        snr2 = (A(1,1)^2*( g(5,401:T+400)*ctranspose(g(5,401:T+400)) + g(6,401:T+400)*ctranspose(g(6,401:T+400)) + g(7,401:T+400)*ctranspose(g(7,401:T+400)) + g(8,401:T+400)*ctranspose(g(8,401:T+400)) ) )/(4*T*2*V);
        SNR2(ppp,nv) = 10*log10(real(snr2));
        snr3 = (A(1,1)^2*( g(9,401:T+400)*ctranspose(g(9,401:T+400)) + g(10,401:T+400)*ctranspose(g(10,401:T+400)) + g(11,401:T+400)*ctranspose(g(11,401:T+400)) + g(12,401:T+400)*ctranspose(g(12,401:T+400)) ) )/(4*T*2*V);
        SNR3(ppp,nv) = 10*log10(real(snr3));

        Bfinal=zeros(K*2,T);
        xh = zeros(4,K);
        Ph = zeros(4,4,K);
        M=4;                                         %no of modes
        sym = [-1 -1; -1 1; 1 -1; 1 1];
        Tr  = [0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25];
        modePr0 = ones(M,K)/M;
        modePr(:,:,1) = modePr0;                       %posterior symbol prob. at time t-1 
        i   = sqrt(-1);

        for t=1:T,
            B   = zeros(K*4, K*4);
            xp1 = zeros(4, 1);
            Pp1 = zeros(4, 4);
            xp2 = zeros(0,1);
            xp3 = zeros(4,1);
            Pp2 = zeros(4, 4);
            P   = zeros(0, 0);
            V1  = zeros(4,4);
            V2  = zeros(4,4);
            H1  = zeros(0, 0);
            mask = [1:1:5];
            if ismember(mod(t,15),mask)
                bhat1 = diag([kron(b(:,t),ones(2,1))]);
                bhat2 = blkdiag(-b(2,t)*eye(2) , b(1,t)*eye(2) , -b(4,t)*eye(2) , b(3,t)*eye(2) , -b(6,t)*eye(2) , b(5,t)*eye(2));
                yt1 = bhat1 * y1(:,t);
                yt2 = bhat2 * y2(:,t);
                Ct1(:,:) = bhat1 * LT * A * bhat1;
                Ct2(:,:) = bhat2 * LT1 * A * bhat2; 
                for k=1:K,
                    xp1 = xh(:,k);
                    Pp1 = Ph(:,:,k) + Q;
                    
                    H = Ct1((k-1)*4+1:4*k,(k-1)*4+1:4*k);
                    if k ~= 1
                        H1= Ct1((k-1)*4+1:4*k,1:(k-1)*4);
                        xp3 = (H1 * xp2);
                        V1 = H1*P*ctranspose(H1);
                    end;
                    V2 = 0;
                    [xh(:,k),Ph(:,:,k),Rlike(1)]=kalman(yt1((k-1)*4+1:k*4),H,xp1,xp3,V,V2,V1,Pp1);
                    xp2((k-1)*4+1:k*4,1) = xh(:,k);
                    P((k-1)*4+1:k*4,(k-1)*4+1:k*4) = Ph(:,:,k);
                end; 
                Rchmse(1:4,(t-1)*2+1) = Rchmse(1:4,(t-1)*2+1) + conj(abs(real(xp2(1:4,1))) - abs(real(g(1:4,(t-1)*2+1+400)))) .* (abs(real(xp2(1:4,1))) - abs(real(g(1:4,(t-1)*2+1+400))));
                Ichmse(1:4,(t-1)*2+1) = Ichmse(1:4,(t-1)*2+1) + conj(abs(imag(xp2(1:4,1))) - abs(imag(g(1:4,(t-1)*2+1+400)))) .* (abs(imag(xp2(1:4,1))) - abs(imag(g(1:4,(t-1)*2+1+400))));
                ghat(1:4,(t-1)*2+1) = xp2(1:4,1);
                xp2 = zeros(0,1);
                P   = zeros(0, 0);
                for k=1:K,
                    xp1 = xh(:,k);
                    Pp1 = Ph(:,:,k) + Q;
                    
                    H = Ct2((k-1)*4+1:4*k,(k-1)*4+1:4*k);
                    if k ~= 1
                        H1= Ct2((k-1)*4+1:4*k,1:(k-1)*4);
                        xp3 = (H1 * xp2);
                        V1 = H1*P*ctranspose(H1);
                    end;
                    V2 = 0;
                    [xh(:,k),Ph(:,:,k),Rlike(1)]=kalman(yt2((k-1)*4+1:k*4),H,xp1,xp3,V,V2,V1,Pp1);
                    xp2((k-1)*4+1:k*4,1) = xh(:,k);
                    P((k-1)*4+1:k*4,(k-1)*4+1:k*4) = Ph(:,:,k);
                end; 
                Rchmse(1:4,t*2) = Rchmse(1:4,t*2) + conj(abs(real(xp2(1:4,1))) - abs(real(g(1:4,t*2+400)))) .* (abs(real(xp2(1:4,1))) - abs(real(g(1:4,t*2+400))));
                Ichmse(1:4,t*2) = Ichmse(1:4,t*2) + conj(abs(imag(xp2(1:4,1))) - abs(imag(g(1:4,t*2+400)))) .* (abs(imag(xp2(1:4,1))) - abs(imag(g(1:4,t*2+400))));
                ghat(1:4,t*2) = xp2(1:4,1);
                Bfinal(:,t)=b(:,t);

                if ismember(mod(t,15),5)
                    xh1 = xh;               %channel c(k,l) for paths p for user k   
                    Ph1 = Ph;             % diagonal terms should not be complex ?!!!!!!!!!!!!                    
                    xh2 = xh;        
                    Ph2 = Ph;                   
                    xh3 = xh;               %channel c(k,l) for paths p for user k   
                    Ph3 = Ph;
                    xh4 = xh;        
                    Ph4 = Ph;
                    [x,y,m1]=intersect(b(1:2,t)',sym,'rows');
                    modePr(:,1,1) = zeros(4,1);
                    modePr(m1,1,1) = 1;
                    [x,y,m1]=intersect(b(3:4,t)',sym,'rows');
                    modePr(:,2,1) = zeros(4,1);
                    modePr(m1,2,1) = 1;
                    [x,y,m1]=intersect(b(5:6,t)',sym,'rows');
                    modePr(:,3,1) = zeros(4,1);
                    modePr(m1,3,1) = 1;

                end;
            else
                for k=1:K, 
                    %FOR TIME i-1
                
                    Tr  = [0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25];
                    sym = [-1 -1; -1 1; 1 -1; 1 1];
                    %calculation of mixing probabilities
                    for i=1:M,
                        for j=1:M,
                            mu(i,j) = Tr(i,j)*modePr(i,k,1);            %prior probabilities
                        end;
                    end;
                    summu=sum(mu);
                    for j=1:M,
                        mu1(:,j)= mu(:,j)/summu(j);
                    end;
                    mu=mu1;
                    
                    if (isnan(mu)~=0)
                        mu
                        pause;
                    end;
                    %MIXING
                    
                    for j=1:M,
                        xh0m(:,j) = xh1(:,k)*mu(1,j) + xh2(:,k)*mu(2,j) + xh3(:,k)*mu(3,j) + xh4(:,k)*mu(4,j);
                        Ph0m(:,:,j)= mu(1,j)*(Ph1(:,:,k) + (xh1(:,k) - xh0m(:,j))*ctranspose(xh1(:,k)-xh0m(:,j)) ) + mu(2,j)*(Ph2(:,:,k) + (xh2(:,k)-xh0m(:,j))*ctranspose(xh2(:,k)-xh0m(:,j))) + mu(3,j)*(Ph3(:,:,k) + (xh3(:,k) - xh0m(:,j))*ctranspose(xh3(:,k)-xh0m(:,j)) ) + mu(4,j)*(Ph4(:,:,k) + (xh4(:,k)-xh0m(:,j))*ctranspose(xh4(:,k)-xh0m(:,j)));
                    end;
                    
                    
                    
                    %Mode matched filtering
                    
                    %Kalman Filter for user k matched to symbols -1 -1
                    j = 1;
                    B((k-1)*4+1:4*k,(k-1)*4+1:4*k) = blkdiag(sym(j,1)*eye(2), sym(j,2)*eye(2));
                    xp1 = xh0m(:,j);   
                    Pp1 = Ph0m(:,:,j) + Q;
                    C(:,:,j) = LT * A * B;     
                    H = C((k-1)*4+1:4*k,(k-1)*4+1:4*k,j);
                    if k ~= 1
                    H1= C((k-1)*4+1:4*k,1:(k-1)*4,j);
                    xp3 = (H1 * xp2);
                    V1 = H1*P*ctranspose(H1);
                    end;
                    
                    [xh1(:,k),Ph1(:,:,k),Rlike(j)]=kalman(y1((k-1)*4+1:k*4,t),H,xp1,xp3,V,V2,V1,Pp1);
                    
                    
                    %Kalman Filter for user k matched to symbol 1
                    j = 2;
                    B((k-1)*4+1:4*k,(k-1)*4+1:4*k) = blkdiag(sym(j,1)*eye(2), sym(j,2)*eye(2));
                    xp1 = xh0m(:,j);   
                    Pp1 = Ph0m(:,:,j) + Q;
                    C(:,:,j) = LT * A * B;     
                    H = C((k-1)*4+1:4*k,(k-1)*4+1:4*k,j);
                    if k ~= 1
                    H1= C((k-1)*4+1:4*k,1:(k-1)*4,j);
                    xp3 = (H1 * xp2);
                    V1 = H1*P*ctranspose(H1);
                    end;
                    [xh2(:,k),Ph2(:,:,k),Rlike(j)]=kalman(y1((k-1)*4+1:k*4,t),H,xp1,xp3,V,V2,V1,Pp1);
                    
                    
                    j=3;
                    B((k-1)*4+1:4*k,(k-1)*4+1:4*k) = blkdiag(sym(j,1)*eye(2), sym(j,2)*eye(2));
                    xp1 = xh0m(:,j);   
                    Pp1 = Ph0m(:,:,j) + Q;
                    C(:,:,j) = LT * A * B;     
                    H = C((k-1)*4+1:4*k,(k-1)*4+1:4*k,j);
                    if k ~= 1
                    H1= C((k-1)*4+1:4*k,1:(k-1)*4,j);
                    xp3 = (H1 * xp2);
                    V1 = H1*P*ctranspose(H1);
                    end;
                    [xh3(:,k),Ph3(:,:,k),Rlike(j)]=kalman(y1((k-1)*4+1:k*4,t),H,xp1,xp3,V,V2,V1,Pp1);
                    
                    
                    
                    %Kalman Filter for user k matched to symbol 1
                    j=4;
                    B((k-1)*4+1:4*k,(k-1)*4+1:4*k) = blkdiag(sym(j,1)*eye(2), sym(j,2)*eye(2));
                    xp1 = xh0m(:,j);   
                    Pp1 = Ph0m(:,:,j) + Q;
                    C(:,:,j) = LT * A * B;     
                    H = C((k-1)*4+1:4*k,(k-1)*4+1:4*k,j);
                    if k ~= 1
                    H1= C((k-1)*4+1:4*k,1:(k-1)*4,j);
                    xp3 = (H1 * xp2);
                    V1 = H1*P*ctranspose(H1);
                    end;
                    [xh4(:,k),Ph4(:,:,k),Rlike(j)]=kalman(y1((k-1)*4+1:k*4,t),H,xp1,xp3,V,V2,V1,Pp1);
                    
                    Rlike=real(Rlike);
                    Rsumlike=sum(Rlike);
                    if (Rsumlike~=0),
                        Rlike=Rlike/Rsumlike;
                    end;
                    like=Rlike;
                    if(like==0),
                        like
                        pause;
                    end;
                    %                 like=real(like);
                    
                    %Calculation of mode probabilities
                    for j=1:M,
                        modePr(j,k,2)= like(j) * ( Tr(1,j)* modePr(1,k,1) + Tr(2,j)* modePr(2,k,1) + Tr(3,j)* modePr(3,k,1) + Tr(4,j)* modePr(4,k,1) );     %posterior probabilities
                    end;
                    summode=sum(modePr(:,k,2));
                    modePr(:,k,2)=modePr(:,k,2)/summode;
                    [a,bk] = max(modePr(:,k,2));
                    mode = zeros(4,1); 
                    mode(bk,1)  = 1;
                    if (isnan(modePr(:,k,2))~=0)
                        modePr(:,k,2)
                        pause;
                    end;
                    %hard cancellation
                    temp = [xh1(:,k) xh2(:,k) xh3(:,k) xh4(:,k) ];
                    %hard cancellation
                    B((k-1)*4+1:k*4,(k-1)*4+1:k*4)=blkdiag(sym(bk,1)*eye(2),sym(bk,2)*eye(2));
                    xp2((k-1)*4+1:k*4,1) = temp(:,bk);
                    temp = [Ph1(:,:,k) Ph2(:,:,k) Ph3(:,:,k) Ph4(:,:,k) ];
                    P((k-1)*4*L+1:k*4*L,(k-1)*4*L+1:k*4*L) = temp(:,(bk-1)*4*L+1:bk*4*L);
                end;
                Rchmse(1:4,(t-1)*2+1) = Rchmse(1:4,(t-1)*2+1) + conj(abs(real(xp2(1:4,1))) - abs(real(g(1:4,(t-1)*2+1+400)))) .* (abs(real(xp2(1:4,1))) - abs(real(g(1:4,(t-1)*2+1+400))));
                Ichmse(1:4,(t-1)*2+1) = Ichmse(1:4,(t-1)*2+1) + conj(abs(imag(xp2(1:4,1))) - abs(imag(g(1:4,(t-1)*2+1+400)))) .* (abs(imag(xp2(1:4,1))) - abs(imag(g(1:4,(t-1)*2+1+400))));
                ghat(1:4,(t-1)*2+1) = xp2(1:4,1);
                
                %FOR TIME i
                
                B   = zeros(K*4, K*4);
                xp1 = zeros(4, 1);
                Pp1 = zeros(4, 4);
                xp2 = zeros(0,1);
                xp3 = zeros(4,1);
                Pp2 = zeros(4, 4);
                P   = zeros(0, 0);
                V1  = zeros(4,4);
                V2  = zeros(4,4);
                H1  = zeros(0, 0);
                for k=1:K,
                    
                    Tr  = [0 0 1 0; 1 0 0 0; 0 0 0 1; 0 1 0 0];         %for next time instant all possible transitions are known
                    %                  Tr  = [0.00001 0.00001 0.99997 0.00001; 0.99997 0.00001 0.00001 0.00001; 0.00001 0.00001 0.00001 0.99997; 0.00001 0.99997 0.00001 0.00001];         %for next time instant all possible transitions are known
                    %calculation of mixing probabilities
                    
                    for i=1:M,
                        for j=1:M,
                            mu(i,j) = Tr(i,j)*modePr(i,k,2);              %prior probabilities, mode prob is from the previous time instant in ST matrix
                        end;
                    end;
                    %MIXING
                    for j=1:M,
                        xh0m(:,j) = xh1(:,k)*mu(1,j) + xh2(:,k)*mu(2,j) + xh3(:,k)*mu(3,j) + xh4(:,k)*mu(4,j);
                        Ph0m(:,:,j)= mu(1,j)*(Ph1(:,:,k) + (xh1(:,k) - xh0m(:,j))*ctranspose(xh1(:,k)-xh0m(:,j)) ) + mu(2,j)*(Ph2(:,:,k) + (xh2(:,k)-xh0m(:,j))*ctranspose(xh2(:,k)-xh0m(:,j))) + mu(3,j)*(Ph3(:,:,k) + (xh3(:,k) - xh0m(:,j))*ctranspose(xh3(:,k)-xh0m(:,j)) ) + mu(4,j)*(Ph4(:,:,k) + (xh4(:,k)-xh0m(:,j))*ctranspose(xh4(:,k)-xh0m(:,j)));
                    
                    end;
                    
                    
                    %Mode matched filtering
                    
                    %Kalman Filter for user k matched to symbols -1 -1
                    j = 1;
                    B((k-1)*4+1:4*k,(k-1)*4+1:4*k) = blkdiag(sym(j,1)*eye(2), sym(j,2)*eye(2));
                    xp1 = xh0m(:,j);   
                    Pp1 = Ph0m(:,:,j) + Q;
                    C(:,:,j) = LT1 * A * B;     
                    H = C((k-1)*4+1:4*k,(k-1)*4+1:4*k,j);
                    if k ~= 1
                        H1= C((k-1)*4+1:4*k,1:(k-1)*4,j);
                        xp3 = (H1 * xp2);
                        V1 = H1*P*ctranspose(H1);
                    end;
                    [xh1(:,k),Ph1(:,:,k),Rlike(j)]=kalman(y2((k-1)*4+1:k*4,t),H,xp1,xp3,V,V2,V1,Pp1);
                    
                    
                    
                    %Kalman Filter for user k matched to symbol 1
                    j = 2;
                    B((k-1)*4+1:4*k,(k-1)*4+1:4*k) = blkdiag(sym(j,1)*eye(2), sym(j,2)*eye(2));
                    xp1 = xh0m(:,j);   
                    Pp1 = Ph0m(:,:,j) + Q;
                    C(:,:,j) = LT1 * A * B;     
                    H = C((k-1)*4+1:4*k,(k-1)*4+1:4*k,j);
                    if k ~= 1
                        H1= C((k-1)*4+1:4*k,1:(k-1)*4,j);
                        xp3 = (H1 * xp2);
                        V1 = H1*P*ctranspose(H1);
                    end;
                    [xh2(:,k),Ph2(:,:,k),Rlike(j)]=kalman(y2((k-1)*4+1:k*4,t),H,xp1,xp3,V,V2,V1,Pp1);
                    
                    j=3;
                    B((k-1)*4+1:4*k,(k-1)*4+1:4*k) = blkdiag(sym(j,1)*eye(2), sym(j,2)*eye(2));
                    xp1 = xh0m(:,j);   
                    Pp1 = Ph0m(:,:,j) + Q;
                    C(:,:,j) = LT1 * A * B;     
                    H = C((k-1)*4+1:4*k,(k-1)*4+1:4*k,j);
                    if k ~= 1
                        H1= C((k-1)*4+1:4*k,1:(k-1)*4,j);
                        xp3 = (H1 * xp2);
                        V1 = H1*P*ctranspose(H1);
                    end;
                    [xh3(:,k),Ph3(:,:,k),Rlike(j)]=kalman(y2((k-1)*4+1:k*4,t),H,xp1,xp3,V,V2,V1,Pp1);
                    
                    %Kalman Filter for user k matched to symbol 1
                    j=4;
                    B((k-1)*4+1:4*k,(k-1)*4+1:4*k) = blkdiag(sym(j,1)*eye(2), sym(j,2)*eye(2));
                    xp1 = xh0m(:,j);   
                    Pp1 = Ph0m(:,:,j) + Q;
                    C(:,:,j) = LT1 * A * B;     
                    H = C((k-1)*4+1:4*k,(k-1)*4+1:4*k,j);
                    if k ~= 1
                        H1= C((k-1)*4+1:4*k,1:(k-1)*4,j);
                        xp3 = (H1 * xp2);
                        V1 = H1*P*ctranspose(H1);
                    end;
                    [xh4(:,k),Ph4(:,:,k),Rlike(j)]=kalman(y2((k-1)*4+1:k*4,t),H,xp1,xp3,V,V2,V1,Pp1);
                    
                    Rlike=real(Rlike);
                    Rsumlike=sum(Rlike);
                    if (Rsumlike~=0),
                        Rlike=Rlike/Rsumlike;
                    end;
                    like=Rlike;
                    if(like==0),
                        like
                        pause;
                    end;
                    %                 like=real(like);
                    
                    
                    %Re-arranging state estimates for the M modes
                    Tr1  = [0 0 1 0; 1 0 0 0; 0 0 0 1; 0 1 0 0];         %for next time instant all possible transitions are known
                    
                    mode_states=[xh1(:,k) xh2(:,k) xh3(:,k) xh4(:,k)];
                    for j=1:M,
                        stat(:,Tr1(:,j)'*[1 2 3 4]')=mode_states(:,j);
                    end;
                    xh1(:,k)=stat(:,1);
                    xh2(:,k)=stat(:,2);
                    xh3(:,k)=stat(:,3);
                    xh4(:,k)=stat(:,4);
                    %Re-arranging variance estimates for the M modes
                    mode_var=[Ph1(:,:,k) Ph2(:,:,k) Ph3(:,:,k) Ph4(:,:,k)];
                    for j=1:M,
                        i=Tr1(:,j)'*[1 2 3 4]';
                        mvar(:,(i-1)*4+1:i*4)=mode_var(:,(j-1)*4+1:j*4);
                    end;
                    Ph1(:,:,k)=mvar(:,1:4);
                    Ph2(:,:,k)=mvar(:,5:8);
                    Ph3(:,:,k)=mvar(:,9:12);
                    Ph4(:,:,k)=mvar(:,13:16);
                    
                    %Calculation of mode probabilities and re-arranging 
                    for j=1:M,
                        modePr(Tr1(:,j)'*[1 2 3 4]',k,1)= like(j) * ( Tr(1,j)* modePr(1,k,2) + Tr(2,j)* modePr(2,k,2) + Tr(3,j)* modePr(3,k,2) + Tr(4,j)* modePr(4,k,2) );     %posterior probabilities
                    end;
                    summode=sum(modePr(:,k,1));
                    modePr(:,k,1)=modePr(:,k,1)/summode;
                    [a,bk] = max(modePr(:,k,1));
                    if (isnan(modePr(:,k,1))~=0)
                        modePr(:,k,1)
                        pause;
                    end; 
                    temp = [xh1(:,k) xh2(:,k) xh3(:,k) xh4(:,k) ];
                    %hard cancellation
                    B((k-1)*4+1:k*4,(k-1)*4+1:k*4)=blkdiag(sym(bk,1)*eye(2),sym(bk,2)*eye(2));
                    Bfinal((k-1)*2+1:k*2,t)=sym(bk,:)';
                    xp2((k-1)*4+1:k*4,1) = temp(:,bk);
                    temp = [Ph1(:,:,k) Ph2(:,:,k) Ph3(:,:,k) Ph4(:,:,k) ];
                    P((k-1)*4*L+1:k*4*L,(k-1)*4*L+1:k*4*L) = temp(:,(bk-1)*4*L+1:bk*4*L);
                end;
                Rchmse(1:4,t*2) = Rchmse(1:4,t*2) + conj(abs(real(xp2(1:4,1))) - abs(real(g(1:4,t*2+400)))) .* (abs(real(xp2(1:4,1))) - abs(real(g(1:4,t*2+400))));
                Ichmse(1:4,t*2) = Ichmse(1:4,t*2) + conj(abs(imag(xp2(1:4,1))) - abs(imag(g(1:4,t*2+400)))) .* (abs(imag(xp2(1:4,1))) - abs(imag(g(1:4,t*2+400))));
                ghat(1:4,t*2) = xp2(1:4,1);
                xh(:,1) = xp2(1:4,1);
                xh(:,2) = xp2(5:8,1);
                xh(:,3) = xp2(9:12,1);
                Ph(:,:,1) = P(1:4,1:4);
                Ph(:,:,2) = P(5:8,5:8);
                Ph(:,:,3) = P(9:12,9:12);
                
            end;
        end;
        %T_b=fix(50000/20)*10;
        BER1fusr10(ppp,nv)=sum(abs(b(1,1:T)-Bfinal(1,1:T)))/(2*(T));
        BER1fusr11(ppp,nv)=sum(abs(b(2,1:T)-Bfinal(2,1:T)))/(2*(T));
        BER1fusr20(ppp,nv)=sum(abs(b(3,1:T)-Bfinal(3,1:T)))/(2*(T));
        BER1fusr21(ppp,nv)=sum(abs(b(4,1:T)-Bfinal(4,1:T)))/(2*(T));
        BER1fusr30(ppp,nv)=sum(abs(b(5,1:T)-Bfinal(5,1:T)))/(2*(T));
        BER1fusr31(ppp,nv)=sum(abs(b(6,1:T)-Bfinal(6,1:T)))/(2*(T));%         
    end;
    Rmseqp = Rchmse/(ppp);
    Imseqp = Ichmse/(ppp);
    mseqp = Rmseqp + sqrt(-1)*Imseqp;
    msef(1,nv) = (mseqp(1,:)*ctranspose(mseqp(1,:)))/(g(1,401:2*T+400)*ctranspose(g(1,401:2*T+400)));
    msef(2,nv) = (mseqp(2,:)*ctranspose(mseqp(2,:)))/(g(2,401:2*T+400)*ctranspose(g(2,401:2*T+400)));
    msef(3,nv) = (mseqp(3,:)*ctranspose(mseqp(3,:)))/(g(3,401:2*T+400)*ctranspose(g(3,401:2*T+400)));
    msef(4,nv) = (mseqp(4,:)*ctranspose(mseqp(4,:)))/(g(4,401:2*T+400)*ctranspose(g(4,401:2*T+400)));
    nv;
    
% end;  
usr10=mean(BER1fusr10)    
usr11=mean(BER1fusr11)
usr20=mean(BER1fusr20)    
usr21=mean(BER1fusr21)
usr30=mean(BER1fusr30)    
usr31=mean(BER1fusr31)
snrf=mean(SNR)
snrf2=mean(SNR2)
snrf3=mean(SNR3)
msef(1,:)
msef(2,:)
msef(3,:)
msef(4,:)
ghat(1:4,1:300)
norm1=msef(1,:)
norm2=msef(2,:)
norm3=msef(3,:)
norm4=msef(4,:)
cest=ghat(1:4,1:150)

% if 1
%     figure,subplot(221), 
%     plot([1:100],real(ghat(1,1:100)),'--',[1:100],real(g(1,1+400:100+400)));   % show the pattern of transmitted symbols
%     grid,title('channel tracking');  xlabel('time'),ylabel('g11')
%     
%     subplot(222),
%     plot([1:100],real(ghat(2,1:100)),'--',[1:100],real(g(2,1+400:100+400)));   % show the pattern of transmitted symbols
%     grid,title('channel tracking');  xlabel('time'),ylabel('g12')
%     
%     subplot(223),
%     plot([1:100],real(ghat(3,1:100)),'--',[1:100],real(g(3,1+400:100+400)));   % show the pattern of transmitted symbols
%     grid,title('channel tracking');  xlabel('time'),ylabel('g21')
%     
%     subplot(224),
%     plot([1:100],real(ghat(4,1:100)),'--',[1:100],real(g(4,1+400:100+400)));   % show the pattern of transmitted symbols
%     grid,title('channel tracking');  xlabel('time'),ylabel('g22')
% end
% x=1:999
% figure,plot(x,real(g(1,2:T)),'-r',x,real(chan(1,2:T)),'.');
% figure,plot(real(snr),usr10(:),'-r',real(snr),usr11(:),'.');
%figure,semilogy(real(snrf),usr10(:),'*-k',real(snrf),usr11(:),'o--k',real(snrf2),usr20(:),'x:k',real(snrf2),usr21(:),'s-.k',real(snrf3),usr30(:),'+-k',real(snrf3),usr31(:),'.:k');

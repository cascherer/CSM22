%Code to generate results in
%
%Dissipativity and Integral Quadratic Constraints
%Tailored computational robustness tests for complex interconnections
%Carsten Scherer

%Parametric system
pa=ureal('pa',0);
G=ss([-1 1;0 -1],[0;pa],[-2 1],0);
sys=[1;G]*G*[1 1];
alv=linspace(0.9,1,100);
%%
%mutlipliers
Pbr=[1 0;0 -1];
%Pbr=[0 1;1  0];
%%Pbr=[0 1;1 -2];


%Static multiplier 
ov=[];
for al=alv;
    p.sys=usubs(sys,'pa',al);
    p.type='D';
    p.uc=[1 1];
    p.psi=tf(1);
    p.P0=Pbr;
    s=robinv(p);a
    ov=[ov s.ov];    
end;
ovs=ov;


%Dynamic multiplier 
s=zpk('s');
psi=ss( ((s+.9)/(s+1))^2);

ov=[];
for al=alv;
    p.sys=usubs(sys,'pa',al);
    p.uc=[1 1];
    p.psi=[1;psi];
    p.P0=Pbr;
    s=robinv(p);
    ov=[ov s.ov];    
end;
ovd=ov;

%Dynamic multiplier with contraint in [26]
ov=[];
for al=alv;
    p.sys=usubs(sys,'pa',al);
    p.uc=[1 1];
    p.psi=[1;psi];
    p.P0=Pbr;    
    s=rinvc(p);
    ov=[ov s.ov];    
end;
ovc=ov;

save resFig5
%% %%%%%%%%%%%%%%%%%
load resFig5
figure(1)
clf
ymax=10;
max([ovs./ovd])
[ovc./ovd]
ovd

plot(alv,ovs,'b',alv,ovd,'r',alv,ovc,'k');grid on;
xlabel('Parameter $\alpha$','interpreter','latex');
ylabel('$\sqrt{{\rm trace}(Y)}$','interpreter','latex')

a=axis;
a(1)=0.9;
a(4)=ymax;axis(a)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=robinv(p);
%Robust ellipsoidal invariance
%
%Information about system:
%p.sys:  System with state-space description
%p.uc:   Uncertainty channel dimension: [number of components of z, number of components of w]
%
%Information about multiplier:
%p.type: 'zf' or 'D'
%p.psi:  Column vector of Zames Falb filters as in paper if p.tpye='zf'
%        General stable basis if p.tpye='D'
%p.P0:   2x2 index matrix: Sector bound if p.type='zf'
%        Region of uncertainy if p.type='D'
%
%Output
%s.ov:   Square root of minimial trace of Y
%s.lmi:  Checkset information about LMI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set default type to p.type='D'
if ~isfield(p,'type')
    p.type='nom';
end;
if ~strcmp(p.type,'zf') & ~strcmp(p.type,'D')
    p.type='nom';
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Dimensions of channels
[k,m]=size(p.sys);
nz=p.uc(1);ne=k-nz;
nw=p.uc(2);nd=m-nw;
[kp,mp]=size(p.psi);

%Build system from [w;d] to [z;w;e;d]
%Output dimensions: nz, nw, ne, nd
sys=[p.sys;eye(m)];  %system from [w;d] to [z;e;w;d]
sys=sys([1:nz k+(1:nz) nz+(1:ne) end-nd+1:end],:);
n=size(sys.a,1);


lmi=[];
%%% Nominal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(p.type,'nom');
    Psi=ss([],[],[],[]);
    P=[];
    Ze=zeros(n);
    sys=sys(nz+nw+(1:(nd+ne)),nw+1:end);
end;
%%% Build Zames Falb multiplier %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(p.type,'zf');
    if nw~=nz;
        error('Uncertainty channel must be square.');
    end
    Ppas=[0 1;1 0];
    %conic combination coefficients
    la=sdpvar(1,kp+1);
    
    %Build Psi and P
    %Sector multiplier always included
    Psi=ss([],[],[],[blkdiag(eye(nz),eye(nz))]);
    P=la(1)*p.P0;
    lmi=lmi+[la(1)>0];
    
    for j=1:kp;
        Psi=[Psi;blkdiag(p.psi(j,:)*eye(nz),eye(nz))];
        P=blkdiag(P,la(j+1)*Ppas);
        lmi=lmi+[la(j+1)>0];
    end
    %Terminal cost zero
    Ze=zeros(size(Psi.a,1)+n);
end;

%%% Build D multiplier %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(p.type,'D');
    if nw~=nz | mp~=nz;
        error('Uncertainty and multipliers do not fit.');
    end
    Psi=blkdiag(p.psi,p.psi);
    
    %Constraint on M
    [Ap,Bp,Cp,Dp]=ssdata(p.psi);
    np=size(Ap,1);
    Xp=sdpvar(np,np);
    M=sdpvar(kp,kp);
    lmi=lmi+[[Ap'*Xp+Xp*Ap Xp*Bp;Bp'*Xp zeros(mp)]-[Cp Dp]'*M*[Cp Dp]<0];
    
    %Multiplier parameter
    P=kron(p.P0,M);
    %Terminal cost matrix
    Z=kron(p.P0,Xp);
    Ze=blkdiag(Z,zeros(n));
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Build common ingredients
%Psi-filtered system
%Output dimensions: nv, ne, nd
nv=size(Psi,1);
sysf=blkdiag(Psi,eye(ne+nd))*sys;

%Extract information about filtered system
[Af,Bf,Cf,Df]=ssdata(sysf);
nf=size(Af,1);
[kf,mf]=size(sysf);
Cfe=Cf(nv+(1:ne),:);
if norm(Df(nv+(1:ne),:))>0;error('[De Ded] should vanish.');end;

%Special performance index matrix
Pp=[zeros(ne) zeros(ne,nd);zeros(nd,ne) -eye(nd)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define LMI system
Xf=sdpvar(nf,nf);
Y=sdpvar(ne,ne);

Pe=blkdiag(P,Pp);

%Build main LMI
lmi=lmi+[[Af'*Xf+Xf*Af Xf*Bf;Bf'*Xf zeros(mf)]+[Cf Df]'*Pe*[Cf Df]<0];
%Coupling and performance LMI
lmi=lmi+[[Y Cfe;Cfe' Xf+Ze]>0];

%Solve LMI: Call Matlab LMI sover from Yalmip
%           A direct implementation is more efficient
h=solvesdp(lmi,trace(Y),sdpsettings('solver','lmilab'));
if h.problem>0;
    s.ov=inf;
else
    s.lmi=checkset(lmi);
    s.ov=sqrt(trace(double(Y)));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=rinvc(p);
%Robust ellipsoidal invariance based on [26]
%
%Information about system:
%p.sys:  System with state-space description
%p.uc:   Uncertainty channel dimension: [number of components of z, number of components of w]
%
%Information about multiplier:
%p.psi:  Filter with state-space description for multiplier 
%p.P0:   2x2 index matrix for uncertainty description 

%Dimensions of channels
[k,m]=size(p.sys);
nz=p.uc(1);ne=k-nz;
nw=p.uc(2);nd=m-nw;

%Build system from [w;d] to [z;w;e;d]
%Output dimensions: nz, nw, ne, nd
sys=[p.sys;eye(m)];  %system from [w;d] to [z;e;w;d]
sys=sys([1:nz k+(1:nz) nz+(1:ne) end-nd+1:end],:);
n=size(sys.a,1);

%Build Psi from psi
Psi=blkdiag(p.psi,p.psi);
nv=size(Psi,1);

%Build Psi-filtered system
%Output dimensions: nv, ne, nd
sysf=blkdiag(Psi,eye(ne+nd))*sys;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extra information about filtered system
[Af,Bf,Cf,Df]=ssdata(sysf);
nf=size(Af,1);
[kf,mf]=size(sysf);
Cfe=Cf(nv+(1:ne),:);
if norm(Df(nv+(1:ne),:))>0;error('[De Ded] should vanish.');end;

%Special performance index matrix
Pp=[zeros(ne) zeros(ne,nd);zeros(nd,ne) -eye(nd)];

%Extra information about psi
[Ap,Bp,Cp,Dp]=ssdata(p.psi);
np=size(Ap,1);
[kp,mp]=size(p.psi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define LMI variables
Xf=sdpvar(nf,nf);
Xp=sdpvar(np,np);
M=sdpvar(kp,kp);
Y=sdpvar(ne,ne);

%Multiplier parameter
P=kron(p.P0,M);
Pe=blkdiag(P,Pp);

%Terminal cost matrix
Z=kron(p.P0,Xp);
Ze=blkdiag(Z,zeros(n));

%Build and solve LMI
lmi=[];
%Constraint on multiplier
lmi=lmi+[[Ap'*Xp+Xp*Ap Xp*Bp;Bp'*Xp zeros(mp)]+[Cp Dp]'*M*[Cp Dp]>0];
%Main LMI 
lmi=lmi+[[Af'*Xf+Xf*Af Xf*Bf;Bf'*Xf zeros(mf)]+[Cf Df]'*Pe*[Cf Df]<0];
%Coupling and performance LMI
lmi=lmi+[[Y Cfe;Cfe' Xf-0*Ze]>0];
lmi=lmi+[M>0];
%Solve LMI 
h=solvesdp(lmi,trace(Y),sdpsettings('solver','lmilab'));
if h.problem>0;
    s.ov=inf;
else
    s.lmi=checkset(lmi);
    s.ov=sqrt(trace(double(Y)));    
end;
end
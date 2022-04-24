%Code to generate results in
%
%Dissipativity and Integral Quadratic Constraints
%Tailored computational robustness tests for complex interconnections
%Carsten Scherer


%parametric system
pa=ureal('pa',0);
G0=ss([-3 -2;1 0],[1;0],[0 -1],0);
G=G0*pa;
sys=[1;1]*G*[1 1];

alv=linspace(5,50,450);
pole1=10;
pole2=100;
%%
%mutlipliers
%Pbr=[1 0;0 -1];
%Ppr=[0 1;1  0];
Psec=[0 1;1 -2];


%Nominal invariance
ov=[];
for al=alv;
    p.sys=usubs(sys,'pa',al);
    p.type='nom';
    p.uc=[1 1];
    p.psi=[];
    p.P0=Psec;
    s=robinv(p);
    ov=[ov s.ov];
end;
ovn=ov;

%Static multiplier
ov=[];
for al=alv;
    p.sys=usubs(sys,'pa',al);
    p.type='zf';
    p.uc=[1 1];
    p.psi=[];
    p.P0=Psec;
    s=robinv(p);
    ov=[ov s.ov];
end;
ovs=ov;

%Dynamic multiplier with pole1
s=zpk('s');
psi1=ss(-pole1,pole1,-1,1);
ov=[];
for al=alv;
    p.sys=usubs(sys,'pa',al);
    p.type='zf';
    p.uc=[1 1];
    p.psi=[psi1];
    p.P0=Psec;
    s=robinv(p);
    ov=[ov s.ov];
end;
ovd1=ov;

%Dynamic multiplier with pole2
psi2=ss(-pole2,pole2,-1,1);
ov=[];
for al=alv;
    p.sys=usubs(sys,'pa',al);
    p.type='zf';
    p.uc=[1 1];
    p.psi=[psi2];
    p.P0=Psec;
    s=robinv(p);
    ov=[ov s.ov];
end;
ovd2=ov;

save resFig4
%% %%%%%%%%%%%%%%%%%
load resFig4
Fontsize=12;
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesFontSize',Fontsize);
figure(1)
clf
ymax=20;
plot(alv,ovn,'k',alv,ovs,'b',alv,ovd1,'r',alv,ovd2,'g');grid on;
p0=alv(min(find(ovs==Inf)));
if ~isempty(p0);
    h=line([p0;p0],[0;ymax],'Color','b','LineStyle',':');
end;
p1=alv(min(find(ovd1==Inf)));
if ~isempty(p1);
    h=line([p1;p1],[0;ymax],'Color','r','LineStyle',':');
end;
p2=alv(min(find(ovd2==Inf)));
if ~isempty(p2);
    h=line([p2;p2],[0;ymax],'Color','g','LineStyle',':');
end;

xlabel('Parameter $\alpha$','interpreter','latex');
ylabel('$\sqrt{{\rm trace}(Y)}$','interpreter','latex')
a=axis;a(4)=ymax;axis(a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
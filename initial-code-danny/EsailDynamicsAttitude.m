clear all
clc
% parpool(4);
    
%% CONSTANTS
muearth=398600.4e9; %gravitational constant of Earth
R0=6378.14e3; %radius of Earth's surface
musun=132712439935.5e9; %gravitation constant of Sun
Rearth=149.5979e9; %mean distance from Sun to Earth
Rmars=227.9423e9; %mean distance from Sun to Mars
thetaplot=0:2*pi/10000:2*pi*9999/10000; %variable theta values for polar plot
Rearthplot=Rearth*ones(1,length(thetaplot)); %Earth distance line for polar plot
Rmarsplot=Rmars*ones(1,length(thetaplot)); %Mars distance line for polar plot
Vearth=sqrt(musun/Rearth); %velocity of Earth in circular orbit around Sun
Vesc=sqrt(2*muearth/R0); %escape velocity of vehicle launched from R0 on Earth
ms=287; %mass of the main spacecraft body
omegax=0; %constant spin rate of spacecraft
KP=0.1; %proportional gain
KI=2; %integral gain
KD=0.5; %derivative gain
R1s=5; %distance from CM to tether 
M=8; %number of tethers
N=10; %number of nodes per tether
L=1000; %unstretched length of tethers (m)
Llink=L/(N-1); %length of tether links (unstretched distance between nodes) (m)
E=72e9; %modulus of elasticity of tethers (N/m^2)
rhoVm=2700; %volumetric mass density of tethers (kg/m^3)
rhoLm=1e-5; %linear mass density of tethers (kg/m)
Aeff=rhoLm/rhoVm; %effective cross-section area of tethers (m^2)
kl=E*Aeff/Llink; %linear stiffness of tethers (N/m)
nul=3e-6; %linear damping ratio of tethers
cl=nul*sqrt(E*Aeff*rhoLm); %linear damping coefficient of tethers (N-s/m)
I=Aeff^2/4/pi; %area moment of inertia of tethers (m^4)
kb=E*I/Llink; %bending stiffness of tethers (N-m/rad)
nub=1e-4; %bending damping ratio of tethers
cb=nub*sqrt(E*Aeff^2*rhoLm*Llink^2/12/pi); %bending damping coefficient of tethers (N-m-s/rad)
massend=rhoLm*L; %end mass on tethers (kg)
Ixx=6000; %x-axis mass moment of inertia of main spacecraft body (N-m-s^2/rad)
It=3000; %tangential (y- and z- axis) mass moment of inertia of main spacecraft body (N-m-s^2/rad)
IB=[Ixx,0,0;0,It,0;0,0,It]; %mass moment of inertia tensor of main spacecraft body (N-m-s^2/rad)
epsilon0=8.85418782e-12; %permittivity of free space constant (F/m)
k=1/(4*pi*epsilon0); %electrostatic constant (m/F)
e=1.60217662e-19; %charge of single electron and proton (C)
kTe=8*e; %thermal energy of solar wind plasma (J)
n0earth=7.3e6; %number density of solar wind plasma at Earth orbit (m^-3)
risetime=20e-3; %time from zero to steady-state voltage on tethers (s)
lambdaDeearth=sqrt(epsilon0*kTe/n0earth/e^2); %Debye length at Earth orbit (m)
V0ss=5.6e3; %steady-state voltage on tethers (V)
rw=1e-3; %effective tether wire radius (m)
rhoLqss=2*pi*epsilon0*V0ss/log(Llink/rw); %linear charge density of tethers at steady-state voltage (C/m)
q=rhoLqss*Llink; %standard point charge at tether nodes (C)
EMerror=1e-6; %maximum magnitude of EM force between two point charges of charge q to be neglected (N)
% syms x
% EMscalefun=@(x) (EMerror*lambdaDeearth^2/k/q^2*x^2-(1+x)*exp(-x)); 
% EMscale=fzero(EMscalefun,3.2132); %characteristic length of EM grid cells for approximate EM force calculations (Debye lengths (m))
EMscale=3.2132;
K=3.09; %solar wind force constant
mp=1.6726219e-27; %mass of a proton (kg)
vSW=4e5; %constant velocity of the solar wind plasma relative to Sun-centered inertial frame (m/s)
maxstep=1e6; %maximum number of allowed simulated time steps

%% VARIABLE PRE-ALLOCATION

%vectors of length 3 in matrices represent directions with the first ...
%representing the x direction, second the y, and third the z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL PRE-ALLOCATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmnsO=cell(M,N); %position vector from main spacecraft body CM to tether nodes in inertial (O) coordinates (m)
rmnsB=cell(M,N); %position vector from main spacecraft body CM to tether nodes in body (B) coordinates (m)
rmnsBold=cell(M,N); %previous step postion vector of tether nodes
rmnp1nB=cell(M,N-1); %position vector from one node to the next on the same tether in B coordinates (m/s)
dx=cell(M,N-1); %length of stretch in tether link
Ovmnp1nB=cell(M,N-1); %velocity of node relative to previous node on same tether; Sun-centered inertial (O) frame, B coordinates (m/s)
rhatmnp1nB=cell(M,N-1); %unit position vector from one node to the next on the same tether; B coordinates
tetherpos=cell(1,M); %rmnsB arranged in 1XM cell of 3*N matrices for 3-D plot (m)
BvmnsB=cell(M,N); %velocity of tether node relative to CM; B frame and coordinates (m/s)
OvmnsO=cell(M,N); %velocity of tether node relative to CM; O frame, O coordinates (m/s)
OvmnsB=cell(M,N); %velocity of tether node relative to CM; O frame, B coordinates (m/s)
BamnsB=cell(M,N); %acceleration of tether node relative to CM; B frame and coordinates (m/s^2)
OamnoB=cell(M,N); %acceleration of tether node relative to Sun; O frame, B coordinates (m/s^2)
OamnsO=cell(M,N); %acceleration of tether node relative to CM; O frame, O coordinates (m/s^2)
OamnsB=cell(M,N); %acceleration of tether node relative to CM; O frame, B coordinates (m/s^2)
massmn=cell(M,N); %mass of node (kg)
qmn=cell(M,N); %charge of node (C)
FGmn=cell(M,N); %gravitational force on node (N)
FTmn=cell(M,N); %tension force on node (N)
FTlink=cell(M,N-1); %tension force on tether link (N)
FTlinkmag=cell(M,N-1); %magnitude of tension force on tether link (N)
TBmn=cell(M,N-2); %bending torque on node (N/m)
phimn=cell(M,N-2); %bending angle of link relative to previous link (rad)
phimndot=cell(M,N-2); %bending angle rotation rate of link relative to previous link (rad/s)
FBmn=cell(M,N); %bending force on node (N)
FSWmn=cell(M,N); %solar wind force on node (N)
FSWlink=cell(M,N-1); %solar wind force on link (N)
FEMmn=cell(M,N); %electrostatic force 
FSUMmn=cell(M,N); %sum of forces on node (N)
EMnodeindex=cell(M,N); %index of EM cell location of each tether node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix Pre-Allocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmnsO(:,:)={zeros(3,1)};
rmnsB(:,:)={zeros(3,1)};
rmnsBold(:,:)={zeros(3,1)};
rmnp1nB(:,:)={zeros(3,1)};
dx(:,:)={zeros(3,1)};
Ovmnp1nB(:,:)={zeros(3,1)};
rhatmnp1nB(:,:)={zeros(3,1)};
tetherpos(:)={zeros(3,N)};
BvmnsB(:,:)={zeros(3,1)};
OvmnsO(:,:)={zeros(3,1)};
OvmnsB(:,:)={zeros(3,1)};
BamnsB(:,:)={zeros(3,1)};
OamnoB(:,:)={zeros(3,1)};
OamnsO(:,:)={zeros(3,1)};
OamnsB(:,:)={zeros(3,1)};
FGmn(:,:)={zeros(3,1)};
FTmn(:,:)={zeros(3,1)};
FTlink(:,:)={zeros(3,1)};
FTlinkmag(:,:)={0};
TBmn(:,:)={zeros(3,1)};
phimn(:,:)={zeros(3,1)};
phimndot(:,:)={zeros(3,1)};
FBmn(:,:)={zeros(3,1)};
FSWmn(:,:)={zeros(3,1)};
FSWlink(:,:)={zeros(3,1)};
FEMmn(:,:)={zeros(3,1)};
FSUMmn(:,:)={zeros(3,1)};
EMnodeindex(:,:)={zeros(3,1)};
rsoOvec=zeros(3,maxstep+1); %time history of position of spacecraft relative to Sun (m)
thetazorbitvec=zeros(1,maxstep+1); %time history of anomaly of spacecraft relative to Sun (rad)
thetazorbitavgerrvec=zeros(1,maxstep+1); %time history of average relative error of anomaly
Rsovec=zeros(1,maxstep+1); %time history of distance of spacecraft from Sun (m)
Rsoavgerrvec=zeros(1,maxstep+1); %time history of average relative error of distance of spacecraft
Rsoavgvec=zeros(1,maxstep+1); %time history of average distance of spacecraft from Sun (m)
thetavec=zeros(3,maxstep+1); %time history of Euler angle values (rad)
thetadotvec=zeros(3,maxstep+1); %time history of Euler angle velocity values (rad/s)
thetaavgerrvec=zeros(3,maxstep+1); %time history of average absolute error of Euler angles from target angles (rad)
thetadotavgerrvec=zeros(3,maxstep+1); %time history of average absolute error of Euler angle rates (rad/s)
thetaztargetvec=zeros(1,maxstep+1); %time history of target value of yaw (z-axis Euler angle) (rad)
thetazdottargetvec=zeros(1,maxstep+1); %time history of of target value of yaw rate (rad/s)
tauctrlvec=zeros(3,maxstep+1); %time history of control torque on spacecraft body (N-m)
tauctrlavgvec=zeros(3,maxstep+1); %time history of average control torque on spacecraft body (N-m)
f=zeros(3,M*N+2,53); %state variables with M*N+1 being spacecraft CM position, ...
    %M*N+2 being euler angles, and the third vector being the order, ...
    %starting from first order at 1, then second order (first derivative) at 2, etc.
fp1=zeros(3,M*N+2,53); %state variables at the next time step
fp1star=zeros(3,M*N+2,53); %predicted state variables at the next time step (previous iteration)
dtp=zeros(1,M*N+2); %maximun acceptable time step for desired accuracy of each state variable
t=zeros(1,maxstep+1); %elapsed simulation time at each time step
dt=zeros(1,maxstep); %time step size at each time step
rm1sB=zeros(3,M); %positoin of 1st tether nodes relative to CM; B coordinates (m)
FSUMm1=zeros(3,M); %sum of the forces on the 1st tether node, excluding that due to spacecraft body
f1error=zeros(1,M*N+2); %absolute change in 0th derivative of state variables at current iteration relative to previous iteration
f2error=zeros(1,M*N+2); %absolute change in 1st derivative of state variables at current iteration relative to previous iteration
f3error=zeros(1,M*N+2); %absolute change in 2nd derivative of state variables at current iteration relative to previous iteration
itervec=zeros(1,maxstep); %time step history of number of iterations til convergence
cputimevec=zeros(1,maxstep); %time step history of total CPU runtime of each step (s)
cputimesumvec=zeros(1,maxstep); %time step history of total elapsed CPU time from beginning of time step loop to the end of each step (s)
iteravgvec=zeros(1,maxstep); %time step history of average number of iterations til convergence
cputimeavgvec=zeros(1,maxstep); %time step history of average CPU time of each step (s)
dtavgvec=zeros(1,maxstep); %time step history of average step size
thrust=zeros(3,maxstep+1); %time history of thrust on spacecraft body (N)
thrustavg=zeros(3,maxstep+1); %time history of average thrust on spacecraft body (N)
thrustmag=zeros(1,maxstep+1); %time history of magnitude of thrust (N) 
thrustavgmag=zeros(1,maxstep+1); %time history of average magnitude of thrust (N)
velRMS=zeros(1,maxstep+1); %time history of RMS of velocity of tether nodes in B frame
accRMS=zeros(1,maxstep+1); %time history of RMS of acceleration of tether nodes in B frame
SEth=zeros(1,M*N+1); %theoretical total energy of tether nodes and spacecraft body
SEact=zeros(1,M*N+1); %actual computed total energy of tether nodes and spacecraft body
kc=zeros(1,M*N+1); %stiffness of artificial correcting spring for energy conservation
cc=zeros(1,M*N+1); %damping coefficient of artificial correcting spring
xc=zeros(1,M*N+1); %equilibrium position of artificial correcting spring
accelck=zeros(3,M*N+1); %acceleration due to artificial spring stiffness
trueanomaly=zeros(1,maxstep+1); %time history of theoretical true anomaly for comparison to computational
r=zeros(1,maxstep+1); %time history of theoretical distance from sun for comparison to computational
R1Nserrvec=zeros(1,maxstep+1);
R1Nsavgerrvec=zeros(1,maxstep+1);
maxFTvec=zeros(1,maxstep+1);
maxFTavgvec=zeros(1,maxstep+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIAL CONDITIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORBITAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsoO=[Rearth;0;0]; %position of spacecraft relative to Sun; O coordinates (m)
OvsoO=[0;0;0]; %velocity of spacecraft relative to Sun; O frame, O coordinates (m/s)
% a=musun*norm(rsoO)/(2*musun-norm(OvsoO)^2*norm(rsoO)); %semi-major axis of transfer ellipse with no additional propulsion (m)
% eccentricity=1-norm(rsoO)/a; %eccentricity of transfer ellipse with no additional propulsion
% Me=sqrt(musun/a^3); %mean motion of transfer ellipse with no additional propulsion (rad/s)
% r(1)=norm(rsoO);
% trueanomaly(1)=0;
SEth(M*N+1)=1/2*norm(OvsoO)^2-musun/norm(rsoO);
SEact(M*N+1)=SEth(M*N+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATTITUDE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetazorbit=atan(rsoO(2)/rsoO(1)); %measured counterclockwise around z-axis in O frame (rad)
thetazdotorbit=1/norm(rsoO)*dot(cross(rsoO/norm(rsoO),OvsoO),[0;0;1]); %tangential velocity over radial distance in O frame (rad/s)
thetaytarget=0; %controller target Euler angle in y direction (B frame) (rad)
thetaydottarget=0; %controller target Euler angle velocity in y direction (B frame) (rad/s)
thetaztarget=thetazorbit; %controller target Euler angle in z direction (B frame) (rad)
thetazdottarget=thetazdotorbit; %controller target Euler angle velocity in z direction (B frame) (rad/s)
theta=[0;thetaytarget;thetaztarget]; % Euler angle vector (rad)
thetadot=[omegax;thetaydottarget;thetazdottarget]; %Euler angle velocity vector (rad/s)
omegaOB=[1,0,-sin(theta(2));0,cos(theta(1)),cos(theta(2))*sin(theta(1)); ... %angular velocity of body w.r.t. inertial frame (rad/s)
    0,-sin(theta(1)),cos(theta(1))*cos(theta(2))]*thetadot; 
taub=[0;0;0]; %external torque on spacecraft body (N-m)
omegaOBcross=[0,-omegaOB(3),omegaOB(2);omegaOB(3),0,-omegaOB(1); ... %angular velocity cross product matrix (rad/s)
    -omegaOB(2),omegaOB(1),0];
tauctrl=[dot(omegaOBcross*IB*omegaOB-taub,[1;0;0]);-KP*(theta(2)- ... %control torque on spacecraft body (N-m)
    thetaytarget)-KD*(thetadot(2)-thetaydottarget);-KP*(theta(3)- ...
    thetaztarget)-KD*(thetadot(3)-thetazdottarget)];
alphaOB=IB\(tauctrl-omegaOBcross*IB*omegaOB); %angular acceleration of body w.r.t. inertial frame (rad/s^2)
thetaddot=[0,cos(theta(1))*tan(theta(2))*thetadot(1)+sin(theta(1)) ... %Euler angle acceleration (rad/s^2)
    *thetadot(2)/cos(theta(2))^2,-sin(theta(1))*tan(theta(2)) ...
    *thetadot(1)+cos(theta(1))*thetadot(2)/cos(theta(2))^2;0, ...
    -sin(theta(1))*thetadot(1),-cos(theta(1))*thetadot(1); 0, ...
    cos(theta(1))*thetadot(1)/cos(theta(2))+sin(theta(1)) ...
    *tan(theta(2))*thetadot(2),(-sin(theta(1))*thetadot(1) ...
    +cos(theta(1))*tan(theta(2))*thetadot(2))/cos(theta(2))] ...
    * omegaOB + [1,sin(theta(1))*tan(theta(2)),cos(theta(1)) ...
    *tan(theta(2));0,cos(theta(1)),-sin(theta(1));0,sin(theta(1)) ...
    /cos(theta(2)),cos(theta(1))/cos(theta(2))]*alphaOB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COORDINATE TRANSFORMATION MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cBO=[cos(theta(3))*cos(theta(2)),cos(theta(2))*sin(theta(3)),-sin(theta(2)) ... %O to B frame coordinate transformation matrix
    ;cos(theta(3))*sin(theta(1))*sin(theta(2))-cos(theta(1))*sin(theta(3)), ...
    cos(theta(1))*cos(theta(3))+sin(theta(1))*sin(theta(2))*sin(theta(3)), ...
    cos(theta(2))*sin(theta(1));cos(theta(1))*cos(theta(3))*sin(theta(2))+ ...
    sin(theta(1))*sin(theta(3)),-cos(theta(3))*sin(theta(1))+cos(theta(1))* ...
    sin(theta(2))*sin(theta(3)),cos(theta(1))*cos(theta(2))];
cOB=cBO.'; %B to O frame coordinate transformation matrix
cBS=[cos(theta(3)-thetazorbit)*cos(theta(2)),cos(theta(2))*sin(theta(3)- ... %S to B frame coordinate transformation matrix
    thetazorbit),-sin(theta(2));cos(theta(3)-thetazorbit)*sin(theta(1))* ...
    sin(theta(2))-cos(theta(1))*sin(theta(3)-thetazorbit),cos(theta(1))* ...
    cos(theta(3)-thetazorbit)+sin(theta(1))*sin(theta(2))*sin(theta(3)- ...
    thetazorbit),cos(theta(2))*sin(theta(1));cos(theta(1))*cos(theta(3)- ...
    thetazorbit)*sin(theta(2))+sin(theta(1))*sin(theta(3)-thetazorbit), ...
    -cos(theta(3)-thetazorbit)*sin(theta(1))+cos(theta(1))*sin(theta(2))* ...
    sin(theta(3)-thetazorbit),cos(theta(1))*cos(theta(2))];
cSB=cBS.'; %B to S frame coordinate transformation matrix
cSO=[cos(thetazorbit),sin(thetazorbit),0;-sin(thetazorbit), ... %O to S frame coordinate transformation matrix
    cos(thetazorbit),0;0,0,1];
cOS=cSO.'; %S to O frame coordinate transformation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TETHER AND SPACECRAFT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:M
    massmn{m,N}=1/2*rhoLm*Llink+massend; %last node consists of end mass and half of last tether segment
    for n=2:N-1
        massmn{m,n}=rhoLm*Llink; %sum of half of both surrounding tether segments
    end
    massmn{m,1}=1/2*rhoLm*Llink; %half of first tether segment
end

%gravity
for m=1:M
    for n=1:N
        %always in -x direction in S frame; approximated by position of
        %spacecraft CM (CM to tether node distances negligible)
        FGmn{m,n}=-massmn{m,n}*musun/norm(rsoO)^2*(cSB*[1;0;0]);
    end
end
FGs=-ms*musun/norm(rsoO)^2*(cSB*[1;0;0]); %main spacecraft body gravitational force

%thrust
thrust(:,1)=zeros(3,1);
thrustavg(:,1)=thrust(:,1);
thrustmag(1)=norm(thrust(:,1));
thrustavgmag(1)=norm(thrustavg(:,1));

%spacecraft body acceleration
OasoB=(FGs+thrust(:,1))/ms;
OasoO=cBO*OasoB;

%tether node kinematics
for m=1:M %position
    %first nodes equally spaced radially around spacecraft CM in Y-Z plane
    rmnsB{m,1}=[0;R1s*cos(2*pi*(m-1)/M);R1s*sin(2*pi*(m-1)/M)];
    for n=2:N
            % unstretched length plus length stretched by centripetal force in
            % straight line with CM and 1st node
            rmnsB{m,n}=rmnsB{m,n-1}+(Llink+(massmn{m,n}*omegaOB(1)^2* ...
                (R1s+Llink))/(kl-massmn{m,n}*omegaOB(1)^2))*rmnsB{m,1}/R1s;           
    end
end
R1Ns=norm(rmnsB{1,N});
for m=1:M %velocity in O frame
    for n=1:N
        OvmnsB{m,n}=BvmnsB{m,n}+cross(omegaOB,rmnsB{m,n});
    end
end
for m=1:M %position relative to previous node
    for n=1:N-1
        rmnp1nB{m,n}=rmnsB{m,n+1}-rmnsB{m,n};
    end
end
for m=1:M %stretch in tether link
    for n=1:N-1
        dx{m,n}=norm(rmnp1nB{m,n})-Llink;
    end
end
for m=1:M %linear tension forces on links
    for n=1:N-1
        if dx{m,n}>0
            FTlink{m,n}=-(kl*dx{m,n}+cl*dot(Ovmnp1nB{m,n},rhatmnp1nB{m,n}))* ...
                rhatmnp1nB{m,n};
        else
            FTlink{m,n}=zeros(3,1);
        end
        FTlinkmag{m,n}=norm(FTlink{m,n});
    end
end
maxFTvec(1)=max(max(cell2mat(FTlinkmag)));
maxFTavgvec(1)=maxFTvec(1);
for m=1:M %acceleration in O frame
    for n=1:N
        OamnsB{m,n}=BamnsB{m,n}+cross(alphaOB,rmnsB{m,n})+2*cross( ...
            omegaOB,BvmnsB{m,n})+cross(omegaOB,cross(omegaOB,rmnsB{m,n}));
        OamnoB{m,n}=OamnsB{m,n}+OasoB;
    end
end
rmnsBold=rmnsB;
rmnsO=rmnsB;
OvmnsO=OvmnsB;
OamnsO=OamnsB;

%EM grid parameters
lambdaDe=lambdaDeearth; %Debye length
EMcellsize=EMscale*lambdaDe; %side length of cubic EM grid cell (m)
EMgridsize=ceil(max(max(max(abs(cell2mat(rmnsB)))))/EMcellsize+2)*2; %side length of cubic EM grid in # of cells
EMcellindex=cell(EMgridsize,EMgridsize,EMgridsize); %index of tether nodes in each EM cell
for m=1:M %indexing node content in EM cells and vice-versa
    for n=1:N
        for i=1:3
            EMnodeindex{m,n}(i)=EMgridsize/2+ceil(rmnsB{m,n}(i)/EMcellsize);
        end
        EMcellindex{EMnodeindex{m,n}(1),EMnodeindex{m,n}(2),EMnodeindex{m,n}(3)} ...
            =[EMcellindex{EMnodeindex{m,n}(1),EMnodeindex{m,n}(2),EMnodeindex{m,n}(3)} ...
            ;[m,n]];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:M %tether nodes (second index of 1:M*N)
    for n=1:N
        f(:,(m-1)*N+n,1)=rmnsO{m,n};
        f(:,(m-1)*N+n,2)=OvmnsO{m,n};
        f(:,(m-1)*N+n,3)=OamnsO{m,n};
    end
end
%main spacecraft body position (second index of M*N+1)
f(:,M*N+1,1)=rsoO;
f(:,M*N+1,2)=OvsoO;
f(:,M*N+1,3)=OasoO;
%attitude (Euler angles) (second index of M*N+2)
f(:,M*N+2,1)=theta;
f(:,M*N+2,2)=thetadot;
f(:,M*N+2,3)=thetaddot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vectors for plotting
rsoOvec(:,1)=rsoO;
thetazorbitvec(1)=thetazorbit;
thetazorbitavgerrvec(1)=0;
Rsovec(1)=norm(rsoO);
Rsoavgvec(1)=Rsovec(1);
Rsoavgerrvec(1)=0;
thetavec(:,1)=theta;
thetadotvec(:,1)=thetadot;
thetaztargetvec(1)=thetaztarget;
thetazdottargetvec(1)=thetazdottarget;
thetaavgerrvec(:,1)=thetavec(:,1)-[0;thetaytarget;thetaztarget];
thetadotavgerrvec(:,1)=thetadotvec(:,1)-[0;thetaydottarget;thetazdottarget];
tauctrlvec(:,1)=tauctrl;
tauctrlavgvec(:,1)=tauctrlvec(:,1);
tauctrlRMS=sqrt(mean(tauctrlvec(:,1).'.^2)).';
thetayRMS=sqrt(mean((thetavec(2,1)-thetaytarget).^2));
thetazRMS=sqrt(mean((thetavec(3,1)-thetaztargetvec(1)).^2));
thetaydotRMS=sqrt(mean((thetadotvec(2,1)-thetaydottarget).^2));
thetazdotRMS=sqrt(mean((thetadotvec(3,1)-thetazdottargetvec(1)).^2));
t(1)=0;
for m=1:M
    for n=1:N
        tetherpos{m}(:,n)=[rmnsB{m,n}(1);rmnsB{m,n}(2);rmnsB{m,n}(3)];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lim=[-size(EMcellindex,1)/2,size(EMcellindex,1)/2]; %axis limits on 2-D tether node plots
tick=-size(EMcellindex,1)/2:1:size(EMcellindex,1)/2; %tick marks on 2-D tether node plots
ticklabel=cell(1,length(tick)); %tick labels on 2-D tether node plots
for i=1:floor(length(tick)/5)
    ticklabel(1+(i-1)*5)={tick(1+(i-1)*5)};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-D TETHER NODE PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure(1);
plot3(tetherpos{1}(1,:),tetherpos{1}(2,:),tetherpos{1}(3,:),'Color','r', ...
    'LineStyle','-','Marker','.','markersize',6)
for m=2:M
    hold on
    plot3(tetherpos{m}(1,:),tetherpos{m}(2,:),tetherpos{m}(3,:),'Color','b', ...
        'LineStyle','-','Marker','.','markersize',6)
end
title '3-D Tether Plot'
xlabel 'x (m)'
ylabel 'y (m)'
zlabel 'z (m)'
axishandle=gca;
set(axishandle,'xlim',[-EMgridsize*EMcellsize/2,EMgridsize*EMcellsize/2], ...
    'ylim',[-EMgridsize*EMcellsize/2,EMgridsize*EMcellsize/2], ...
    'zlim',[-EMgridsize*EMcellsize/2,EMgridsize*EMcellsize/2]);
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Y-Z PLANE TETHER NODE PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig2=figure(2);
clf
plot(tetherpos{1}(2,:)/EMcellsize,tetherpos{1}(3,:)/EMcellsize,'Color','r', ...
    'LineStyle','-','Marker','.','markersize',6)
for m=2:M
    hold on
    plot(tetherpos{m}(2,:)/EMcellsize,tetherpos{m}(3,:)/EMcellsize,'Color', ...
        'b','LineStyle','-','Marker','.','markersize',6)
end
title 'Y-Z Plane Tether Plot'
xlabel 'y (Cell Lengths)'
ylabel 'z (Cell Lengths)'
axishandle=gca;
set(axishandle,'xlim',lim,'ylim',lim,'xtick',tick,'ytick',tick,'xgrid','on', ...
    'ygrid','on','xticklabel',ticklabel,'yticklabel',ticklabel)
drawnow
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%X-Z PLANE TETHER NODE PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig3=figure(3);
clf
plot(tetherpos{1}(1,:)/EMcellsize,tetherpos{1}(3,:)/EMcellsize,'Color','r', ...
    'LineStyle','-','Marker','.','markersize',6)
for m=2:M
    hold on
    plot(tetherpos{m}(1,:)/EMcellsize,tetherpos{m}(3,:)/EMcellsize,'Color', ...
        'b','LineStyle','-','Marker','.','markersize',6)
end
title 'X-Z Plane Tether Plot'
xlabel 'x (Cell Lenghts)'
ylabel 'z (Cell Lengths)'
axishandle=gca;
set(axishandle,'xlim',lim,'ylim',lim,'xtick',tick,'ytick',tick,'xgrid','on', ...
    'ygrid','on','xticklabel',ticklabel,'yticklabel',ticklabel)
drawnow
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORBIT PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig4=figure(4);
clf
polarplot(thetazorbitvec(1),Rsovec(1),'linestyle','-','marker','.','markersize',6)
hold on
polarplot(thetaplot,Rearthplot,'linestyle','--','color','k')
hold on
polarplot(thetaplot,Rmarsplot,'linestyle','--','color','k')
title 'Orbit'
rlim([0,Rmars*1.1])
drawnow
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE ALLOCATION FOR OTHER PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%theta plot
fig5=figure(5);

%tauctrl plot
fig6=figure(6);

%thrust plot
fig7=figure(7);

%time plot
fig8=figure(8);

%B frame RMS plot
fig9=figure(9);

%elapsed time
fig10=figure(10);

%stability plot/max tension plot
fig11=figure(11);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TIME STEP LOOP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME STEP LOOP CONSTANTS AND PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dttol=1e7; %error tolerance for calculating time step size
sstol=5e-4; %error tolerance for steady-state convergence
ss=[0,0]; %steady-state vector with 1st value representing steady-state
%convergence condition (0=false, 1=true) and 2nd value being the time step
%at which it occurs
rise=1; %risetime condition (0=false, 1=true)
rev=0; %number of completed orbit revolutions (thetazorbit)
orbit=0; %index of whether spacecraft is completing an orbit (vs. small deviations in path)
half=1; %index of which half (upper or lower) of orbit spacecraft is in (1=upper, 0=lower)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME STEP FOR LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for step=1:maxstep
        
    starttime=cputime; %start of computational time tracking
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP SIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:M*N+2
        if step==1 %first step size
            dtp(i)=1e-4;
        else
            %recommended step size based on individual state variable
            %calculated by comparing magnitude of 2nd order term to
            %predicted magnitude of 3rd order term (f3p) while limiting
            %the allowed increase in step size from the previous time step
            dtp(i)=100/(5/(dttol*norm(f(:,i,3))/norm(f3p(:,i)))+95/dt(step-1));
            %limiting step size based on stability criterion
            %(CFL=sqrt(E/rhoVm)<=1)
            if dtp(i)>0.99*Llink/sqrt(E/rhoVm)
                dtp(i)=0.99*Llink/sqrt(E/rhoVm);
            end
        end
    end
    dt(step)=min(dtp); %step size based on smallest recommended step size of all state variables
    t(step+1)=t(step)+dt(step); %time value of next step
    
    %forcing time step to end at the end of the transient state of the wire
    %voltage instead of stepping over the start of the steady-state point
    if rise==1 && t(step+1)>risetime
        dt(step)=risetime-t(step);
        t(step+1)=risetime;
        rise=0;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSIENT VOLTAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    if t(step+1)<=risetime
        V0=V0ss*(t(step)+dt(step))/risetime; %voltage rises linearly
        rhoLq=2*pi*epsilon0*V0/log(Llink/rw); %linear chage density
        for m=1:M
            qmn{m,N}=1/2*rhoLq*Llink;
            for n=2:N-1
                qmn{m,n}=rhoLq*Llink;
            end
            qmn{m,1}=1/2*rhoLq*Llink;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMERICAL TIME INTEGRATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTEGRATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:3 %setup for 1st iteration (0th iteration state variable values)
        for j=1:M*N+2
            %previous iteration, next time step value of state variables
            fp1star(i,j,1)=f(i,j,1)+f(i,j,2)*dt(step);
            fp1star(i,j,2)=f(i,j,2);
        end
    end
    
    for iter=1:50 %iteration loop for predictor-corrector method (max of 50 iterations allowed)
        
        for i=1:3 %current iteration state variable values (Taylor series)
            for j=1:M*N+2
                %second term is iter+1 order Taylor series term
                fp1(i,j,1)=fp1star(i,j,1)+1/factorial(iter+1)* ...
                    f(i,j,iter+2)*dt(step)^(iter+1);
                fp1(i,j,2)=fp1star(i,j,2)+1/factorial(iter)*f(i,j,iter+2)* ...
                    dt(step)^(iter);
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATED VARIABLES OF INTEREST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for m=1:M 
            for n=1:N
                rmnsO{m,n}=fp1(:,(m-1)*N+n,1);
                OvmnsO{m,n}=fp1(:,(m-1)*N+n,2);
            end
        end
        rsoO=fp1(:,M*N+1,1);
        OvsoO=fp1(:,M*N+1,2);
        theta=fp1(:,M*N+2,1);
        thetadot=fp1(:,M*N+2,2);
        omegaOB=[1,0,-sin(theta(2));0,cos(theta(1)),cos(theta(2))*sin(theta(1)); ...
         0,-sin(theta(1)),cos(theta(1))*cos(theta(2))]*thetadot;
        thetazorbit=atan(rsoO(2)/rsoO(1));
        thetazdotorbit=1/norm(rsoO)*dot(cross(rsoO/norm(rsoO),OvsoO),[0;0;1]);
        thetaztarget=0;
        thetaztargetvec(step+1)=thetaztarget;
        thetazdottarget=0;
        thetazdottargetvec(step+1)=thetazdottarget;
        % coordinate transformation matrices
        cBO=[cos(theta(3))*cos(theta(2)),cos(theta(2))*sin(theta(3)),-sin(theta(2)) ...
            ;cos(theta(3))*sin(theta(1))*sin(theta(2))-cos(theta(1))*sin(theta(3)), ...
            cos(theta(1))*cos(theta(3))+sin(theta(1))*sin(theta(2))*sin(theta(3)), ...
            cos(theta(2))*sin(theta(1));cos(theta(1))*cos(theta(3))*sin(theta(2))+ ...
            sin(theta(1))*sin(theta(3)),-cos(theta(3))*sin(theta(1))+cos(theta(1))* ...
            sin(theta(2))*sin(theta(3)),cos(theta(1))*cos(theta(2))];
        cOB=cBO.';
        cBS=[cos(theta(3)-thetazorbit)*cos(theta(2)),cos(theta(2))*sin(theta(3)- ...
            thetazorbit),-sin(theta(2));cos(theta(3)-thetazorbit)*sin(theta(1))* ...
            sin(theta(2))-cos(theta(1))*sin(theta(3)-thetazorbit),cos(theta(1))* ...
            cos(theta(3)-thetazorbit)+sin(theta(1))*sin(theta(2))*sin(theta(3)- ...
            thetazorbit),cos(theta(2))*sin(theta(1));cos(theta(1))*cos(theta(3)- ...
            thetazorbit)*sin(theta(2))+sin(theta(1))*sin(theta(3)-thetazorbit), ...
            -cos(theta(3)-thetazorbit)*sin(theta(1))+cos(theta(1))*sin(theta(2))* ...
            sin(theta(3)-thetazorbit),cos(theta(1))*cos(theta(2))];
        cSB=cBS.';
        cSO=[cos(thetazorbit),sin(thetazorbit),0;-sin(thetazorbit), ...
            cos(thetazorbit),0;0,0,1];
        cOS=cSO.';
        
        for m=1:M
            for n=1:N
                rmnsB{m,n}=cBO*rmnsO{m,n};
                OvmnsB{m,n}=cBO*OvmnsO{m,n};
                BvmnsB{m,n}=OvmnsB{m,n}-cross(omegaOB,rmnsB{m,n});
            end
        end        
        for m=1:M
            for n=1:N-1
                rmnp1nB{m,n}=rmnsB{m,n+1}-rmnsB{m,n};
                Ovmnp1nB{m,n}=OvmnsB{m,n+1}-OvmnsB{m,n};
                rhatmnp1nB{m,n}=rmnp1nB{m,n}/norm(rmnp1nB{m,n});
            end
        end
        for m=1:M
            for n=1:N-1
                dx{m,n}=norm(rmnp1nB{m,n})-Llink;
            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% TETHER NODE FORCE CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %parallel computing shit (don't worry about it)
%         spmd
%         switch labindex
            
%             case 1
        for m=1:M %gravitational forces on nodes
            for n=1:N
                FGmn{m,n}=-massmn{m,n}*musun/norm(rsoO)^2*(cBS*[1;0;0]);
                FGmn{m,n}=zeros(3,1);
            end
        end
        FGs=-ms*musun/norm(rsoO)^2*(cBS*[1;0;0]); %gravitational force on spacecraft body
        FGs=zeros(3,1);
        for m=1:M %linear tension forces on links
            for n=1:N-1
                if dx{m,n}>0
                    FTlink{m,n}=-(kl*dx{m,n}+cl*dot(Ovmnp1nB{m,n},rhatmnp1nB{m,n}))* ...
                        rhatmnp1nB{m,n};
                else
                    FTlink{m,n}=zeros(3,1);
                end
                FTlinkmag{m,n}=norm(FTlink{m,n});
            end
        end
        for m=1:M %linear tension forces on nodes
            FTmn{m,N}=FTlink{m,N-1};
            for n=2:N-1
                FTmn{m,n}=FTlink{m,n-1}-FTlink{m,n};
            end
            FTmn{m,1}=-FTlink{m,1};
        end
        
%             case 2
        for m=1:M %bending torques on nodes
            for n=1:N-2
                if norm(cross(rhatmnp1nB{m,n},rhatmnp1nB{m,n+1})) ~= 0
                    phimn{m,n}=acos(dot(rmnp1nB{m,n+1},rhatmnp1nB{m,n})/ ...
                        norm(rmnp1nB{m,n+1}))*cross(rhatmnp1nB{m,n}, ...
                        rhatmnp1nB{m,n+1})/norm(cross(rhatmnp1nB{m,n}, ...
                        rhatmnp1nB{m,n+1}));
                    if isreal(phimn{m,n})==0
                        phimn(m,n)={zeros(3,1)};
                    end
                else
                    phimn(m,n)={zeros(3,1)};
                end
                if norm(cross(rhatmnp1nB{m,n},rhatmnp1nB{m,n+1})) ~= 0
                    Ovmnp1nBplane=Ovmnp1nB{m,n+1}-dot(Ovmnp1nB{m,n+1},cross( ...
                        rhatmnp1nB{m,n},rhatmnp1nB{m,n+1})/norm(cross( ...
                        rhatmnp1nB{m,n},rhatmnp1nB{m,n+1})));
                else
                    Ovmnp1nBplane=Ovmnp1nB{m,n+1};
                end
                Ovmnp1nBperp=Ovmnp1nBplane-dot(Ovmnp1nBplane,rhatmnp1nB{m,n+1});
                phimnp1ndot=cross(rhatmnp1nB{m,n+1},Ovmnp1nBperp)/norm(rmnp1nB{m,n+1});
                if norm(cross(rhatmnp1nB{m,n},rhatmnp1nB{m,n+1})) ~= 0
                    Ovmnm1nBplane=-Ovmnp1nB{m,n}-dot(-Ovmnp1nB{m,n},cross( ...
                        rhatmnp1nB{m,n},rhatmnp1nB{m,n+1})/norm(cross( ...
                        rhatmnp1nB{m,n},rhatmnp1nB{m,n+1})));
                else
                    Ovmn1nBplane=-Ovmnp1nB{m,n};
                end
                Ovmnm1nBperp=Ovmnp1nBplane-dot(Ovmnp1nBplane,-rhatmnp1nB{m,n});
                phimnm1ndot=cross(-rhatmnp1nB{m,n},Ovmnm1nBperp)/norm(rmnp1nB{m,n});
                phimndot{m,n}=phimnp1ndot+phimnm1ndot;
                TBmn{m,n}=-(kb*phimn{m,n}+cb*phimndot{m,n});
            end   
        end
        for m=1:M %bending forces on nodes
            FBmn{m,N}=norm(TBmn{m,N-2})/norm(rmnp1nB{m,N-1})*cross( ...
                TBmn{m,N-2},rmnp1nB{m,N-1})/norm(cross(TBmn{m,N-2}, ...
                rmnp1nB{m,N-1}));
            FBmn{m,N-1}=norm(TBmn{m,N-3})/norm(rmnp1nB{m,N-2})*cross( ...
                TBmn{m,N-3},rmnp1nB{m,N-2})/norm(cross(TBmn{m,N-3}, ...
                rmnp1nB{m,N-2}));
            for n=3:N-2
                FBmn{m,n}=norm(TBmn{m,n-2})/norm(rmnp1nB{m,n-1})*cross( ...
                    TBmn{m,n-2},rmnp1nB{m,n-1})/norm(cross(TBmn{m,n-2}, ...
                    rmnp1nB{m,n-1}))-norm(TBmn{m,n})/norm(-rmnp1nB{m,n})* ...
                    cross(TBmn{m,n},-rmnp1nB{m,n})/norm(cross(TBmn{m,n}, ...
                    -rmnp1nB{m,n}));
            end
            FBmn{m,2}=-norm(TBmn{m,2})/norm(-rmnp1nB{m,2})* ...
                    cross(TBmn{m,2},-rmnp1nB{m,2})/norm(cross(TBmn{m,2}, ...
                    -rmnp1nB{m,2}));
            FBmn{m,1}=-norm(TBmn{m,1})/norm(-rmnp1nB{m,1})* ...
                    cross(TBmn{m,1},-rmnp1nB{m,1})/norm(cross(TBmn{m,1}, ...
                    -rmnp1nB{m,1}));
                for n=1:N
                    if isnan(FBmn{m,n})==1
                        FBmn(m,n)={zeros(3,1)};
                    end
                end
        end
        
%             case 3
        %updated SW parameters
        n0=n0earth*(Rearth/norm(rsoO))^2;
        lambdaDe=sqrt(epsilon0*kTe/n0/e^2);
        r0=2*lambdaDe;
        for m=1:M %SW force on links
            for n=1:N-1
                vSWperp=(vSW-dot(OvsoO,rsoO/norm(rsoO)))*norm(cBS*[1;0;0]- ...
                    dot(cBS*[1;0;0],rhatmnp1nB{m,n})*rhatmnp1nB{m,n});
                if isinf(sqrt(exp(mp*vSWperp^2/e/V0*log(r0/rw))-1))~=1
                    FSWlink{m,n}=K*mp*n0*vSWperp^2*r0*Llink/sqrt(exp(mp* ...
                        vSWperp^2/e/V0*log(r0/rw))-1)*(cBS*[norm(rsoO);0;0]-dot(cBS* ...
                        [norm(rsoO);0;0],rhatmnp1nB{m,n}))/(norm(cBS*[norm(rsoO);0;0]-dot(cBS* ...
                        [norm(rsoO);0;0],rhatmnp1nB{m,n})));
                else
                    FSWlink{m,n}=zeros(3,1);
                end
            end
        end
        for m=1:M %SW force on nodes
            FSWmn{m,N}=1/2*FSWlink{m,N-1};
            for n=2:N-1
                FSWmn{m,n}=1/2*(FSWlink{m,n-1}+FSWlink{m,n});
            end
            FSWmn{m,1}=1/2*FSWlink{m,1};
        end
        
%             case 4
        %EM grid parameters
        EMcellsize=EMscale*lambdaDe;
        EMgridsize=ceil(max(max(max(abs(cell2mat(rmnsB)))))/EMcellsize+2)*2;
        EMcellindex=cell(EMgridsize,EMgridsize,EMgridsize);
        for m=1:M %indexing node content in EM cells and vice-versa
            for n=1:N
                for i=1:3
                    EMnodeindex{m,n}(i)=EMgridsize/2+ceil(rmnsB{m,n}(i)/EMcellsize);
                end
                EMcellindex{EMnodeindex{m,n}(1),EMnodeindex{m,n}(2),EMnodeindex{m,n}(3)} ...
                    =[EMcellindex{EMnodeindex{m,n}(1),EMnodeindex{m,n}(2),EMnodeindex{m,n}(3)} ...
                    ;[m,n]];
            end
        end
        for m=1:M %EM forces on node
            for n=1:N
                FEMmn(m,n)={zeros(3,1)};
                for i=EMnodeindex{m,n}(1)-1:EMnodeindex{m,n}(1)+1
                    for j=EMnodeindex{m,n}(2)-1:EMnodeindex{m,n}(2)+1
                        for k=EMnodeindex{m,n}(3)-1:EMnodeindex{m,n}(3)+1
                            for l=1:size(EMcellindex{i,j,k},1)
                                if m==EMcellindex{i,j,k}(l,1) && n==EMcellindex{i,j,k}(l,2)
                                    FEMmn{m,n}=FEMmn{m,n};
                                else
                                FEMmn{m,n}=FEMmn{m,n}+(k*qmn{m,n}*qmn{EMcellindex{i,j,k}(l,1),EMcellindex{ ...
                                    i,j,k}(l,2)}*(1+norm(rmnsB{m,n}-rmnsB{EMcellindex{i,j,k}(l,1), ...
                                    EMcellindex{i,j,k}(l,2)})/lambdaDe)*exp(-norm(rmnsB{m,n}-rmnsB{ ...
                                    EMcellindex{i,j,k}(l,1),EMcellindex{i,j,k}(l,2)})/lambdaDe)/ ...
                                    norm(rmnsB{m,n}-rmnsB{EMcellindex{i,j,k}(l,1),EMcellindex{i,j,k}(l,2)}) ...
                                    ^2)*(rmnsB{m,n}-rmnsB{EMcellindex{i,j,k}(l,1),EMcellindex{i,j,k}(l,2)})/ ...
                                    norm(rmnsB{m,n}-rmnsB{EMcellindex{i,j,k}(l,1),EMcellindex{i,j,k}(l,2)});
                                end
                            end
                        end
                    end
                end
            end
        end

%         end
%         end
%         %more parallel computing shit
%         FGmn=FGmn{1};
%         FGs=FGs{1};
%         FTlinkmag=FTlinkmag{1};
%         FTmn=FTmn{1};
%         FBmn=FBmn{1};
%         FSWmn=FSWmn{1};
%         FEMmn=FEMmn{1};
%         EMgridsize=EMgridsize{1};
%         EMcellsize=EMcellsize{1};

        for m=1:M %sum of forces on tether nodes
            for n=1:N
                FSUMmn{m,n}=FGmn{m,n}+FTmn{m,n}+FBmn{m,n}+FSWmn{m,n}+FEMmn{m,n};
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TETHER/SPACECRAFT BODY INTERACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for m=1:M
            for n=2:N
                OamnoB{m,n}=FSUMmn{m,n}/massmn{m,n};
            end
        end
        
        %calculating tauctrl,x, alphay, and alphaz
        for m=1:M
            rm1sB(:,m)=rmnsB{m,1};
            FSUMm1(:,m)=FSUMmn{m,1};
        end
        m1=massmn{1,1};
        tauctrl(2)=-((KP*(theta(2)-thetaytarget)+KI*((theta(2)- ...
            thetaytarget)+(f(2,M*N+2,1)-thetaytarget))/ ...
            2*dt(step)+KD*(thetadot(2)-thetaydottarget))*cos(theta(1))+ ...
            (KP*(theta(3)-thetaztarget)+KI*((theta(3) ...
            -thetaztargetvec(step+1))+(f(3,M*N+2,1)-thetaztargetvec(step)))/ ...
            2*dt(step)+KD*(thetadot(3)-thetazdottarget))*sin(theta(1))*cos(theta(2)));
        tauctrl(3)=-((KP*(theta(3)-thetaztarget)+KI*((theta(3) ...
            -thetaztargetvec(step+1))+(f(3,M*N+2,1)-thetaztargetvec(step)))/ ...
            2*dt(step)+KD*(thetadot(3)-thetazdottarget))*cos(theta(1))*cos(theta(2))- ...
            (KP*(theta(2)-thetaytarget)+KI*((theta(2)- ...
            thetaytarget)+(f(2,M*N+2,1)-thetaytarget))/ ...
            2*dt(step)+KD*(thetadot(2)-thetaydottarget))*sin(theta(1)));
        A=[1,-IB(1,2)-m1*sum(rm1sB(1,:).*rm1sB(2,:)),-IB(1,3)-m1* ...
            sum(rm1sB(1,:).*rm1sB(3,:));0,-IB(2,2)+m1*sum(rm1sB(1,:).^2+ ...
            rm1sB(3,:).^2),-IB(2,3)-m1*sum(rm1sB(2,:).*rm1sB(3,:));0, ...
            IB(2,3)-m1*sum(rm1sB(2,:).*rm1sB(3,:)),-IB(3,3)+m1*sum( ...
            rm1sB(1,:).^2+rm1sB(2,:).^2)];
        B=[-IB(1,3)*omegaOB(1)*omegaOB(2)-IB(2,3)*omegaOB(2)^2+IB(1,2)* ...
            omegaOB(1)*omegaOB(3)-IB(2,2)*omegaOB(2)*omegaOB(3)+IB(3,3)* ...
            omegaOB(2)*omegaOB(3)-IB(2,3)*omegaOB(3)^2-m1*sum((FSUMm1(3,:).* ...
            rm1sB(2,:)-FSUMm1(2,:).*rm1sB(3,:)).*(-rm1sB(3,:)*omegaOB(2)+ ...
            rm1sB(2,:)*omegaOB(3)).*(rm1sB(1,:)*omegaOB(1)+rm1sB(2,:)* ...
            omegaOB(2)+rm1sB(3,:)*omegaOB(3)));-tauctrl(2)+IB(1,3)*omegaOB(1)^2+ ...
            IB(2,3)*omegaOB(1)*omegaOB(2)+IB(1,1)*omegaOB(1)*omegaOB(3)- ...
            IB(3,3)*omegaOB(1)*omegaOB(3)+IB(1,2)*omegaOB(2)*omegaOB(3)+ ...
            IB(1,3)*omegaOB(3)^2-m1*sum((FSUMm1(3,:).*rm1sB(1,:)-FSUMm1(1,:).* ...
            rm1sB(3,:)).*(-rm1sB(3,:)*omegaOB(1)+rm1sB(1,:)*omegaOB(3)).* ...
            (rm1sB(1,:)*omegaOB(1)+rm1sB(2,:)*omegaOB(2)+rm1sB(3,:)*omegaOB(3))); ...
            -tauctrl(3)-IB(1,2)*omegaOB(1)^2-IB(1,1)*omegaOB(1)*omegaOB(2)+ ...
            IB(2,2)*omegaOB(1)*omegaOB(2)-IB(1,2)*omegaOB(2)^2+IB(2,3)* ...
            omegaOB(1)*omegaOB(3)-IB(1,3)*omegaOB(2)*omegaOB(3)-m1* ...
            sum((FSUMm1(2,:).*rm1sB(1,:)-FSUMm1(1,:).*rm1sB(2,:)).* ...
            (-rm1sB(2,:)*omegaOB(1)+rm1sB(1,:)*omegaOB(2)).*(rm1sB(1,:)* ...
            omegaOB(1)+rm1sB(2,:)*omegaOB(2)+rm1sB(3,:)*omegaOB(3)))];
        C=A\B;
        tauctrl(1)=C(1);
        alphaOB(2)=C(2);
        alphaOB(3)=C(3);
        
        for m=1:M %acceleration of 1st nodes
            OamnsB{m,1}=cross(alphaOB,rm1sB(:,m))+ ...
                cross(omegaOB,cross(omegaOB,rm1sB(:,m)));
            OamnsO{m,1}=cOB*OamnsB{m,1};
        end
        
        thrust(:,step+1)=zeros(3,1);
        for m=1:M %force on spacecraft from tethers
            thrust(:,step+1)=thrust(:,step+1)+(FSUMm1(:,m)-FGmn{m,1}-m1* ...
                OamnsB{m,1})/(m1/ms*M+1);
        end
        OasoB=(thrust(:,step+1)+FGs)/ms;
        OasoO=cOB*OasoB;
        
        for m=1:M 
            OamnoB{m,1}=OamnsB{m,1}+OasoB;
            for n=2:N
                OamnsB{m,n}=OamnoB{m,n}-OasoB;
                BamnsB{m,n}=OamnsB{m,n}-cross(alphaOB,rmnsB{m,n})-cross(omegaOB,BvmnsB{m,n})- ...
                    -cross(omegaOB,cross(omegaOB,rmnsB{m,n}));
                OamnsO{m,n}=cOB*OamnsB{m,n};
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACCELERATION TERMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        thetaddot=[0,cos(theta(1))*tan(theta(2))*thetadot(1)+sin(theta(1)) ...
        *thetadot(2)/cos(theta(2))^2,-sin(theta(1))*tan(theta(2)) ...
        *thetadot(1)+cos(theta(1))*thetadot(2)/cos(theta(2))^2;0, ...
        -sin(theta(1))*thetadot(1),-cos(theta(1))*thetadot(1); 0, ...
        cos(theta(1))*thetadot(1)/cos(theta(2))+sin(theta(1)) ...
        *tan(theta(2))*thetadot(2),(-sin(theta(1))*thetadot(1) ...
        +cos(theta(1))*tan(theta(2))*thetadot(2))/cos(theta(2))] ...
        * omegaOB + [1,sin(theta(1))*tan(theta(2)),cos(theta(1)) ...
        *tan(theta(2));0,cos(theta(1)),-sin(theta(1));0,sin(theta(1)) ...
        /cos(theta(2)),cos(theta(1))/cos(theta(2))]*alphaOB;
        for m=1:M
            for n=1:N
                fp1(:,(m-1)*N+n,3)=OamnsO{m,n};
            end
        end
        fp1(:,M*N+1,3)=OasoO;
        fp1(:,M*N+2,3)=thetaddot;
     
        %energy conservation correction artificial spring term
%         if step>=2
%             fp1(:,M*N+1,3)=fp1(:,M*N+1,3)+accelck(:,M*N+1);
%         end

        %velocity Verlet acceleration
%         accelavg=(fp1(:,M*N+1,3)+f(:,M*N+1,3))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRECTOR TERM CALCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:3
            for j=1:M*N+2
                f(i,j,iter+3)=fp1(i,j,3)*factorial(iter)/dt(step)^(iter);
                for k=1:iter
                    f(i,j,iter+3)=f(i,j,iter+3)-(1/factorial(k-1)*f(i,j,k+2)* ...
                        dt(step)^(k-1))*factorial(iter)/dt(step)^(iter);
                end
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ITERATION CONVERGENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:M*N+2
%                 f1error(i)=norm(fp1(:,i,1)-fp1star(:,i,1));
%                 f2error(i)=norm(fp1(:,i,2)-fp1star(:,i,2));
%                 f3error(i)=norm(fp1(:,i,3)-fp1star(:,i,3));
                f1error(i)=norm(fp1(:,i,1)-fp1star(:,i,1))/norm(fp1(:,i,1));
                f2error(i)=norm(fp1(:,i,2)-fp1star(:,i,2))/norm(fp1(:,i,2));
                f3error(i)=norm(fp1(:,i,3)-fp1star(:,i,3))/norm(fp1(:,i,3));
        end
        if iter>=2
            if max([max(max(f1error)),max(max(f2error)),max(max(f3error))])<=5e-4/maxstep
                break
            end
        end
        fp1star=fp1;
        if step==1
            break
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STABILITY CORRECTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %velocity Verlet
%     fp1(:,:,2)=f(:,:,2)+accelavg*dt(step);

    %Taylor series energy conservation corrections
%     SEth(M*N+1)=SEth(M*N+1)+dot(cOB*thrust(:,step+1)/ms,fp1(:,M*N+1,1)-f(:,M*N+1,1));
%     SEact(M*N+1)=1/2*norm(fp1(:,M*N+1,2))^2-musun/norm(rsoO);
%     if abs(SEact(M*N+1))>abs(SEth(M*N+1))
%         fp1(:,M*N+1,2)=sqrt(norm(fp1(:,M*N+1,2))^2+2*abs(SEact(M*N+1)-SEth(M*N+1)))* ...
%             fp1(:,M*N+1,2)/norm(fp1(:,M*N+1,2));
%     end
%     if abs(SEact(M*N+1))<abs(SEth(M*N+1))
%         fp1(:,M*N+1,2)=sqrt(norm(fp1(:,M*N+1,2))^2-2*abs(SEact(M*N+1)-SEth(M*N+1)))* ...
%             fp1(:,M*N+1,2)/norm(fp1(:,M*N+1,2));
%     end 
%     if abs(SEact(M*N+1))>abs(SEth(M*N+1))
%         xc(M*N+1)=-musun/(1/2*norm(fp1(:,M*N+1,2))^2-SEth(M*N+1))/ ...
%             dot(cross(cross(rsoO,OvsoO),OvsoO)/norm(cross(cross(rsoO, ...
%             OvsoO),OvsoO)),-rsoO/norm(rsoO));
%     else
%         xc(M*N+1)=musun/(1/2*norm(fp1(:,M*N+1,2))^2-SEth(M*N+1))/ ...
%             dot(cross(cross(rsoO,OvsoO),OvsoO)/norm(cross(cross(rsoO, ...
%             OvsoO),OvsoO)),-rsoO/norm(rsoO));
%     end
%     kc(M*N+1)=2*abs(SEact(M*N+1)-SEth(M*N+1))/xc(M*N+1)^2;
%     accelck(:,M*N+1)=kc(M*N+1)*xc(M*N+1)*cross(cross(rsoO,OvsoO),OvsoO)/ ...
%         norm(cross(cross(rsoO,OvsoO),OvsoO))*sign(dot(OvsoO/norm(OvsoO), ...
%         cOS*[sign(xc(M*N+1));0;0]));
%     if xc(M*N+1)==0
%         accelck(:,M*N+1)=zeros(3,1);
%     end
%     fp1(:,M*N+1,3)=fp1(:,M*N+1,3)+accelck(:,M*N+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT VECTOR SHIT (don't worry about it)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    maxFTvec(step+1)=max(max(cell2mat(FTlinkmag)));
    maxFTavgvec(step+1)=(t(step)*maxFTavgvec(step)+(maxFTvec(step+1)+ ...
        maxFTvec(step))/2*dt(step))/t(step+1);
    Rsovec(step+1)=norm(rsoO);
    Rsoavgerrvec(step+1)=(t(step)*Rsoavgerrvec(step)+((Rsovec(step+1)/ ...
        r(step+1)-1)+(Rsovec(step)/r(step)-1))/2*dt(step))/t(step+1);
    Rsoavgvec(step+1)=(t(step)*Rsoavgvec(step)+(Rsovec(step+1)+ ...
        Rsovec(step))/2*dt(step))/t(step+1);
    thetavec(:,step+1)=theta;
    thetadotvec(:,step+1)=thetadot;
    thetaavgerrvec(:,step+1)=(t(step)*thetaavgerrvec(:,step)+((thetavec(:,step+1)- ...
        [0;thetaytarget;thetaztargetvec(step+1)])+(thetavec(:,step)- ...
        [0;thetaytarget;thetaztargetvec(step)]))/2*dt(step))/t(step+1);
    thetadotavgerrvec(:,step+1)=(t(step)*thetadotavgerrvec(:,step)+((thetadotvec(:,step+1)- ...
        [0;thetaydottarget;thetazdottargetvec(step+1)])+(thetadotvec(:,step)- ...
        [0;thetaydottarget;thetazdottargetvec(step)]))/2*dt(step))/t(step+1);
    tauctrlvec(:,step+1)=tauctrl;
    tauctrlavgvec(:,step+1)=(t(step)*tauctrlavgvec(:,step)+(tauctrlvec(:,step+1)+ ...
        tauctrlvec(:,step))/2*dt(step))/t(step+1);
    tauctrlRMS=sqrt(mean(tauctrlvec(:,1:step+1).'.^2)).';
    thetayRMS=sqrt(mean((thetavec(2,1:step+1)-thetaytarget).^2));
    thetazRMS=sqrt(mean((thetavec(3,1:step+1)-thetaztargetvec(1:step+1)).^2));
    thetaydotRMS=sqrt(mean((thetadotvec(2,1:step+1)-thetaydottarget).^2));
    thetazdotRMS=sqrt(mean((thetadotvec(3,1:step+1)-thetazdottargetvec(1:step+1)).^2));
    for m=1:M
        for n=1:N
            velRMS(step+1)=velRMS(step+1)+norm(BvmnsB{m,n})^2;
        end
    end
    velRMS(step+1)=sqrt(velRMS(step+1)/(M*N));
    for m=1:M
        for n=1:N
            accRMS(step+1)=accRMS(step+1)+norm(BamnsB{m,n})^2;
        end
    end
    accRMS(step+1)=sqrt(accRMS(step+1)/(M*N));
    for m=1:M
        for n=1:N
            tetherpos{m}(:,n)=[rmnsB{m,n}(1);rmnsB{m,n}(2);rmnsB{m,n}(3)];
        end
    end
    lim=[-size(EMcellindex,1)/2,size(EMcellindex,1)/2];
    tick=-size(EMcellindex,1)/2:1:size(EMcellindex,1)/2;
    ticklabel=cell(1,length(tick));
    for i=1:floor(length(tick)/5)
        ticklabel(1+(i-1)*5)={tick(1+(i-1)*5)};
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SHIT (magic happens here)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     spmd
%     switch labindex
%         
%         case 1
    %3-D plot
    set(0,'currentfigure',fig1)
    clf
    plot3(tetherpos{1}(1,:),tetherpos{1}(2,:),tetherpos{1}(3,:),'Color','r', ...
        'LineStyle','-','Marker','.','markersize',6)
    for m=2:M
        hold on
        plot3(tetherpos{m}(1,:),tetherpos{m}(2,:),tetherpos{m}(3,:),'Color','b', ...
            'LineStyle','-','Marker','.','markersize',6)
    end
    title '3-D Tether Plot'
    xlabel 'x (m)'
    ylabel 'y (m)'
    zlabel 'z (m)'
    axishandle=gca;
    set(axishandle,'xlim',[-EMgridsize*EMcellsize/2,EMgridsize*EMcellsize/2], ...
        'ylim',[-EMgridsize*EMcellsize/2,EMgridsize*EMcellsize/2], ...
        'zlim',[-EMgridsize*EMcellsize/2,EMgridsize*EMcellsize/2]);
    drawnow

    %Y-Z plot
    set(0,'currentfigure',fig2)
    clf
    plot(tetherpos{1}(2,:)/EMcellsize,tetherpos{1}(3,:)/EMcellsize,'Color','r', ...
        'LineStyle','-','Marker','.','markersize',6)
    for m=2:M
        hold on
        plot(tetherpos{m}(2,:)/EMcellsize,tetherpos{m}(3,:)/EMcellsize,'Color', ...
            'b','LineStyle','-','Marker','.','markersize',6)
    end
    title 'Y-Z Plane Tether Plot'
    xlabel 'y (Cell Lengths)'
    ylabel 'z (Cell Lengths)'
    axishandle=gca;
    set(axishandle,'xlim',lim,'ylim',lim,'xtick',tick,'ytick',tick,'xgrid','on', ...
        'ygrid','on','xticklabel',ticklabel,'yticklabel',ticklabel)
    drawnow

    %X-Z plot
    set(0,'currentfigure',fig3)
    clf
    plot(tetherpos{1}(1,:)/EMcellsize,tetherpos{1}(3,:)/EMcellsize,'Color','r', ...
        'LineStyle','-','Marker','.','markersize',6)
    for m=2:M
        hold on
        plot(tetherpos{m}(1,:)/EMcellsize,tetherpos{m}(3,:)/EMcellsize,'Color', ...
            'b','LineStyle','-','Marker','.','markersize',6)
    end
    title 'X-Z Plane Tether Plot'
    xlabel 'x (Cell Lenghts)'
    ylabel 'z (Cell Lengths)'
    axishandle=gca;
    set(axishandle,'xlim',lim,'ylim',lim,'xtick',tick,'ytick',tick,'xgrid','on', ...
        'ygrid','on','xticklabel',ticklabel,'yticklabel',ticklabel)
    drawnow

    %orbit plot
    set(0,'currentfigure',fig4)
    clf
    polarplot(thetazorbitvec(1:step+1),Rsovec(1:step+1),'linestyle','-','marker','.','markersize',6)
    hold on
    polarplot(thetaplot,Rearthplot,'linestyle','--','color','k')
    hold on
    polarplot(thetaplot,Rmarsplot,'linestyle','--','color','k')
    title 'Orbit'
    drawnow

%         case 2
    %theta plot
    set(0,'currentfigure',fig5)
    clf
    subplot(2,2,1)
    plot(t(1:step+1),thetavec(2,1:step+1),'marker','.','markersize',6)
    hold on
    plot(t(1:step+1),thetaytarget*ones(1,step+1)+thetaavgerrvec(2,1:step+1),'marker','.','markersize',6)
    hold on
    plot(t(1:step+1),thetaytarget*ones(1,step+1))
    hold on
    plot(t(1:step+1),(thetaytarget+thetayRMS)*ones(1,step+1))
    hold on
    plot(t(1:step+1),(thetaytarget-thetayRMS)*ones(1,step+1))
    title '\theta_{y}'
    xlabel 'Time (s)'
    ylabel '\theta_{y} (rad)'
    lgd=legend('\theta_{y}','Mean Error','Target','+RMS','-RMS');
    set(lgd,'color','none')
    if step>20
        axes('position',[.35,.625,.1,.1])
        plot(t(1:step+1),thetavec(2,1:step+1),'--')
        hold on
        plot(t(1:step+1),thetaytarget*ones(1,step+1)+thetaavgerrvec(2,1:step+1),'--')
        hold on
        plot(t(1:step+1),thetaytarget*ones(1,step+1),'--')
        hold on
        plot(t(1:step+1),(thetaytarget+thetayRMS)*ones(1,step+1),'--')
        hold on
        plot(t(1:step+1),(thetaytarget-thetayRMS)*ones(1,step+1),'--')
        set(gca,'color','none')
        minn=min(thetavec(2,step-20:step+1));
        maxx=max(thetavec(2,step-20:step+1));
        if minn==maxx
            minn=minn-1;
            maxx=maxx+1;
        end
        axishandle=gca;
        set(axishandle,'xlim',[t(step-20),t(step+1)],'ylim',[minn-(maxx-minn)*0.1, ...
            maxx+(maxx-minn)*0.1]);
    end
    subplot(2,2,2)
    plot(t(1:step+1),thetavec(3,1:step+1),'marker','.','markersize',6)
    hold on
    plot(t(1:step+1),thetaztargetvec(1:step+1)+thetaavgerrvec(3,1:step+1),'marker','.','markersize',6)
    hold on
    plot(t(1:step+1),thetaztargetvec(1:step+1))
    hold on
    plot(t(1:step+1),(thetaztargetvec(1:step+1)+thetazRMS))
    hold on
    plot(t(1:step+1),(thetaztargetvec(1:step+1)-thetazRMS))
    title '\theta_{z}'
    xlabel 'Time (s)'
    ylabel '\theta_{z} (rad)'
    lgd=legend('\theta_{z}','Mean Error','Target','+RMS','-RMS');
    set(lgd,'color','none')
    if step>20
        axes('position',[.8,.625,.1,.1])
        plot(t(1:step+1),thetavec(3,1:step+1),'--')
        hold on
        plot(t(1:step+1),thetaztargetvec(1:step+1)+thetaavgerrvec(3,1:step+1),'--')
        hold on
        plot(t(1:step+1),thetaztargetvec(1:step+1),'--')
        hold on
        plot(t(1:step+1),(thetaztargetvec(1:step+1)+thetazRMS),'--')
        hold on
        plot(t(1:step+1),(thetaztargetvec(1:step+1)-thetazRMS),'--')
        set(gca,'color','none')
        minn=min(thetavec(3,step-20:step+1));
        maxx=max(thetavec(3,step-20:step+1));
        if minn==maxx
        minn=minn-1;
        maxx=maxx+1;
        end
        axishandle=gca;
        set(axishandle,'xlim',[t(step-20),t(step+1)],'ylim',[minn-(maxx-minn)*0.1, ...
            maxx+(maxx-minn)*0.1]);
    end
    subplot(2,2,3)
    plot(t(1:step+1),thetadotvec(2,1:step+1),'marker','.','markersize',6)
    hold on
    plot(t(1:step+1),thetaydottarget*ones(1,step+1)+thetadotavgerrvec(2,1:step+1),'marker','.','markersize',6)
    hold on
    plot(t(1:step+1),thetaydottarget*ones(1,step+1))
    hold on
    plot(t(1:step+1),(thetaydottarget+thetaydotRMS)*ones(1,step+1))
    hold on
    plot(t(1:step+1),(thetaydottarget-thetaydotRMS)*ones(1,step+1))
    title 'd\theta_{y}/dt'
    xlabel 'Time (s)'
    ylabel 'd\theta_{y}/dt (rad/s)'
    lgd=legend('d\theta_{y}/dt','Mean Error','Target','+RMS','-RMS');
    set(lgd,'color','none')
    if step>20
        axes('position',[.35,.15,.1,.1])
        plot(t(1:step+1),thetadotvec(2,1:step+1),'--')
        hold on
        plot(t(1:step+1),thetaydottarget*ones(1,step+1)+thetadotavgerrvec(2,1:step+1),'--')
        hold on
        plot(t(1:step+1),thetaydottarget*ones(1,step+1),'--')
        hold on
        plot(t(1:step+1),(thetaydottarget+thetaydotRMS)*ones(1,step+1),'--')
        hold on
        plot(t(1:step+1),(thetaydottarget-thetaydotRMS)*ones(1,step+1),'--')
        set(gca,'color','none')
        minn=min(thetadotvec(2,step-20:step+1));
        maxx=max(thetadotvec(2,step-20:step+1));
        if minn==maxx
            minn=minn-1;
            maxx=maxx+1;
        end    
        axishandle=gca;
        set(axishandle,'xlim',[t(step-20),t(step+1)],'ylim',[minn-(maxx-minn)*0.1, ...
            maxx+(maxx-minn)*0.1]);
    end
    subplot(2,2,4)
    plot(t(1:step+1),thetadotvec(3,1:step+1),'marker','.','markersize',6)
    hold on
    plot(t(1:step+1),thetazdottargetvec(1:step+1)+thetadotavgerrvec(3,1:step+1),'marker','.','markersize',6)
    hold on
    plot(t(1:step+1),thetazdottargetvec(1:step+1))
    hold on
    plot(t(1:step+1),(thetazdottargetvec(1:step+1)+thetazdotRMS))
    hold on
    plot(t(1:step+1),(thetazdottargetvec(1:step+1)-thetazdotRMS))
    title 'd\theta_{z}/dt'
    xlabel 'Time (s)'
    ylabel 'd\theta_{z}/dt (rad/s)'
    lgd=legend('d\theta_{z}/dt','Mean Error','Target','+RMS','-RMS');
    set(lgd,'color','none')
    if step>20
        axes('position',[.8,.15,.1,.1])
        plot(t(1:step+1),thetadotvec(3,1:step+1),'--')
        hold on
        plot(t(1:step+1),thetazdottargetvec(1:step+1)+thetadotavgerrvec(3,1:step+1),'--')
        hold on
        plot(t(1:step+1),thetazdottargetvec(1:step+1),'--')
        hold on
        plot(t(1:step+1),(thetazdottargetvec(1:step+1)+thetazdotRMS),'--')
        hold on
        plot(t(1:step+1),(thetazdottargetvec(1:step+1)-thetazdotRMS),'--')
        set(gca,'color','none')
        minn=min(thetadotvec(3,step-20:step+1));
        maxx=max(thetadotvec(3,step-20:step+1));
        if minn==maxx
            minn=minn-1;
            maxx=maxx+1;
        end    
        axishandle=gca;
        set(axishandle,'xlim',[t(step-20),t(step+1)],'ylim',[minn-(maxx-minn)*0.1, ...
            maxx+(maxx-minn)*0.1]);
    end
    drawnow
    
    %tauctrl plot
    set(0,'currentfigure',fig6)
    clf
    subplot(3,1,1)
    plot(t(1:step+1),tauctrlvec(1,1:step+1),'marker','.','markersize',6)
    hold on
    plot(t(1:step+1),tauctrlavgvec(1,1:step+1),'marker','.','markersize',6)
    hold on
    plot(t(1:step+1),tauctrlRMS(1)*ones(1,step+1))
    hold on
    plot(t(1:step+1),-tauctrlRMS(1)*ones(1,step+1))
    title '\tau_{ctrl,x}'
    xlabel 'Time (s)'
    ylabel '\tau_{ctrl,x} (N-m)'
    lgd=legend('\tau_{ctrl,x}','Mean','+RMS','-RMS');
    set(lgd,'color','none')
    subplot(3,1,2)
    plot(t(1:step+1),tauctrlvec(2,1:step+1),'marker','.','markersize',6)
    hold on
    plot(t(1:step+1),tauctrlavgvec(2,1:step+1),'marker','.','markersize',6)
    hold on
    plot(t(1:step+1),tauctrlRMS(2)*ones(1,step+1))
    hold on
    plot(t(1:step+1),-tauctrlRMS(2)*ones(1,step+1))
    title '\tau_{ctrl,y}'
    xlabel 'Time (s)'
    ylabel '\tau_{ctrl,y} (N-m)'
    lgd=legend('\tau_{ctrl,y}','Mean','+RMS','-RMS');
    set(lgd,'color','none')
    subplot(3,1,3)
    plot(t(1:step+1),tauctrlvec(3,1:step+1),'marker','.','markersize',6)
    hold on
    plot(t(1:step+1),tauctrlavgvec(3,1:step+1),'marker','.','markersize',6)
    hold on
    plot(t(1:step+1),tauctrlRMS(3)*ones(1,step+1))
    hold on
    plot(t(1:step+1),-tauctrlRMS(3)*ones(1,step+1))
    title '\tau_{ctrl,z}'
    xlabel 'Time (s)'
    ylabel '\tau_{ctrl,z} (N-m)'
    lgd=legend('\tau_{ctrl,z}','Mean','+RMS','-RMS');
    set(lgd,'color','none')
    drawnow
    
%         case 3
    %thrust plot
    thrustavg(:,step+1)=(t(step)*thrustavg(:,step)+(thrust(:,step+1)+thrust(:,step))/ ...
        2*dt(step))/t(step+1);
    thrustmag(step+1)=norm(thrust(:,step+1));
    thrustavgmag(step+1)=norm(thrustavg(:,step+1));
    set(0,'currentfigure',fig7)
    clf
    plot(t(1:step+1),thrustmag(1:step+1))
    hold on
    plot(t(1:step+1),thrustavgmag(1:step+1))
    title 'Thrust'
    xlabel 'Time (s)'
    ylabel 'Thrust (N)'
    lgd=legend('Thrust','Mean');
    set(lgd,'color','none')
    if step>20
        axes('position',[.65,.15,.25,.25])
        plot(t(1:step+1),thrustmag(1:step+1),'--')
        hold on
        plot(t(1:step+1),thrustavgmag(1:step+1),'--')
        set(gca,'color','none')
        minn=min(thrustmag(step-20:step+1));
        maxx=max(thrustmag(step-20:step+1));
        if minn==maxx
            minn=minn-1;
            maxx=maxx+1;
        end
        axishandle=gca;
        set(axishandle,'xlim',[t(step-20),t(step+1)],'ylim',[minn-(maxx-minn)*0.1, ...
            maxx+(maxx-minn)*0.1]);
    end
    drawnow
    
    %time and iteration plots
    itervec(step)=iter;
    endtime=cputime;
    cputimevec(step)=endtime-starttime;
    if step==1
        iteravgvec(step)=itervec(step);
        cputimeavgvec(step)=cputimevec(step);
        dtavgvec(step)=dt(step);
    else
        iteravgvec(step)=((step-1)*iteravgvec(step-1)+(itervec(step)+ ...
            itervec(step-1))/2)/step;
        cputimeavgvec(step)=((step-1)*cputimeavgvec(step-1)+(cputimevec(step)+ ...
            cputimevec(step-1))/2)/step;
        dtavgvec(step)=((step-1)*dtavgvec(step-1)+(dt(step)+ ...
            dt(step-1))/2)/step;
    end
    set(0,'currentfigure',fig8)
    clf
    subplot(2,2,1)
    plot(1:step,itervec(1:step))
    hold on
    plot(1:step,iteravgvec(1:step))
    title 'Iterations'
    xlabel 'Time Step'
    ylabel 'No. of Iterations'
    lgd=legend('Iterations','Mean');
    set(lgd,'color','none')
    axishandle=gca;
    set(axishandle,'ylim',[0,max(itervec(1:step))+1])
    if step>20
    axes('position',[.35,.625,.1,.1])
    plot(1:step,itervec(1:step),'--')
    hold on
    plot(1:step,iteravgvec(1:step),'--')
    set(gca,'color','none')
    minn=min(iteravgvec(step-20:step));
    maxx=max(iteravgvec(step-20:step));
    if minn==maxx
        minn=minn-1;
        maxx=maxx+1;
    end
    axishandle=gca;
    set(axishandle,'xlim',[step-20,step+1],'ylim',[minn-(maxx-minn)*0.1, ...
        maxx+(maxx-minn)*0.1]);
    end
    subplot(2,2,2)
    plot(1:step,dt(1:step))
    hold on
    plot(1:step,dtavgvec(1:step))
    title 'Step Size'
    xlabel 'Time Step'
    ylabel 'Step Size (s)'
    lgd=legend('Step Size','Mean');
    set(lgd,'color','none')
    if step>20
    axes('position',[.8,.625,.1,.1])
    plot(1:step,dt(1:step),'--')
    hold on
    plot(1:step,dtavgvec(1:step),'--')
    set(gca,'color','none')
    minn=min(dt(step-20:step));
    maxx=max(dt(step-20:step));
    if minn==maxx
        minn=minn-1;
        maxx=maxx+1;
    end
    axishandle=gca;
    set(axishandle,'xlim',[step-20,step+1],'ylim',[minn-(maxx-minn)*0.1, ...
        maxx+(maxx-minn)*0.1]);
    end
    subplot(2,2,3)
    plot(1:step,cputimevec(1:step))
    hold on
    plot(1:step,cputimeavgvec(1:step))
    title 'CPU Time'
    xlabel 'Time Step'
    ylabel 'CPU Time (s)'
    lgd=legend('CPU Time','Mean');
    set(lgd,'color','none')
    if step>20
    axes('position',[.35,.15,.1,.1])
    plot(1:step,cputimevec(1:step),'--')
    hold on
    plot(1:step,cputimeavgvec(1:step),'--')
    set(gca,'color','none')
    minn=min(cputimeavgvec(step-20:step));
    maxx=max(cputimeavgvec(step-20:step));
    if minn==maxx
        minn=minn-1;
        maxx=maxx+1;
    end
    axishandle=gca;
    set(axishandle,'xlim',[step-20,step+1],'ylim',[minn-(maxx-minn)*0.1, ...
        maxx+(maxx-minn)*0.1]);
    end   
    subplot(2,2,4)
    plot(1:step,dt(1:step)./cputimevec(1:step))
    hold on
    plot(1:step,dtavgvec(1:step)./cputimeavgvec(1:step))
    title 'Simulation/CPU Time Ratio'
    xlabel 'Time Step'
    ylabel '(Simulation Time)/(CPU Time)'
    lgd=legend('Sim/CPU Time','Mean');
    set(lgd,'color','none')
    if step>20
    axes('position',[.8,.15,.1,.1])
    plot(1:step,dt(1:step)./cputimevec(1:step),'--')
    hold on
    plot(1:step,dtavgvec(1:step)./cputimeavgvec(1:step),'--')
    set(gca,'color','none')
    minn=min(dtavgvec(step-20:step)./cputimeavgvec(step-20:step));
    maxx=max(dtavgvec(step-20:step)./cputimeavgvec(step-20:step));
    if minn==maxx
        minn=minn-1;
        maxx=maxx+1;
    end
    axishandle=gca;
    set(axishandle,'xlim',[step-20,step+1],'ylim',[minn-(maxx-minn)*0.1, ...
        maxx+(maxx-minn)*0.1]);
    end  
    drawnow
    
%         case 4
    %B frame RMS plot
    set(0,'currentfigure',fig9)
    clf
    subplot(2,1,1)
    plot(t(1:step+1),velRMS(1:step+1))
    hold on
    plot(t(1:step+1),sstol*ones(1,step+1))
    title 'B Frame Velocity RMS'
    xlabel 'Time (s)'
    ylabel 'Velocity RMS (m/s)'
    lgd=legend('Velocity RMS','Convergence');
    set(lgd,'color','none')
    if step>20
        axes('position',[.8,.625,.1,.1])
        plot(t(1:step+1),velRMS(1:step+1),'--')
        hold on
        plot(t(1:step+1),sstol*ones(1,step+1),'--')
        set(gca,'color','none')
        minn=min(velRMS(step-20:step+1));
        maxx=max(velRMS(step-20:step+1));
        if minn==maxx
            minn=minn-1;
            maxx=maxx+1;
        end
        axishandle=gca;
        set(axishandle,'xlim',[t(step-20),t(step+1)],'ylim',[minn-(maxx-minn)*0.1, ...
            maxx+(maxx-minn)*0.1]);
    end
    subplot(2,1,2)
    plot(t(1:step+1),accRMS(1:step+1))
    hold on
    plot(t(1:step+1),sstol*ones(1,step+1))
    title 'B Frame Acceleration RMS'
    xlabel 'Time (s)'
    ylabel 'Acceleration RMS (m/s^{2})'
    lgd=legend('Acceleration RMS','Convergence');
    set(lgd,'color','none')
    if step>20
        axes('position',[.8,.15,.1,.1])
        plot(t(1:step+1),accRMS(1:step+1),'--')
        hold on
        plot(t(1:step+1),sstol*ones(1,step+1),'--')
        set(gca,'color','none')
        minn=min(accRMS(step-20:step+1));
        maxx=max(accRMS(step-20:step+1));
        if minn==maxx
            minn=minn-1;
            maxx=maxx+1;
        end
        axishandle=gca;
        set(axishandle,'xlim',[t(step-20),t(step+1)],'ylim',[minn-(maxx-minn)*0.1, ...
            maxx+(maxx-minn)*0.1]);
    end
    drawnow
    
    %elapsed time plot
    cputimesumvec(step)=sum(cputimevec(1:step));
    set(0,'currentfigure',fig10)
    clf
    subplot(2,1,1)
    plot(1:step,t(2:step+1))
    title 'Elapsed Simulation Time'
    xlabel 'Time Step'
    ylabel 'Simulation Time (s)'
    subplot(2,1,2)
    plot(1:step,cputimesumvec(1:step))
    title 'Elapsed CPU Time'
    xlabel 'Time Step'
    ylabel 'CPU Time (s)'
    drawnow

    %Max Tension Plot
    set(0,'currentfigure',fig11)
    clf
    plot(t(1:step+1),maxFTvec(1:step+1))
    hold on
    plot(t(1:step+1),maxFTavgvec(1:step+1))
    title 'Maximum Tension'
    xlabel 'Time (s)'
    ylabel 'Tension (N)'
    lgd=legend('Tension','Mean');
    set(lgd,'color','none')
    if step>20
        axes('position',[.65,.15,.25,.25])
        plot(t(1:step+1),maxFTvec(1:step+1),'--')
        hold on
        plot(t(1:step+1),maxFTavgvec(1:step+1),'--')
        set(gca,'color','none')
        minn=min(maxFTvec(step-20:step+1));
        maxx=max(maxFTvec(step-20:step+1));
        if minn==maxx
            minn=minn-1;
            maxx=maxx+1;
        end
        axishandle=gca;
        set(axishandle,'xlim',[t(step-20),t(step+1)],'ylim',[minn-(maxx-minn)*0.1, ...
            maxx+(maxx-minn)*0.1]);
    end
    drawnow
    
%     end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEXT TIME STEP PREP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    f3p=zeros(3,M*N+2);
    for i=1:iter-1
        f3p=f3p+1/factorial(i-1)*f(:,:,i+3)*dt(step)^(i-1);
    end
    
    f(:,:,1:3)=fp1(:,:,1:3);
    rmnsBold=rmnsB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME STEP LOOP CONVERGENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rise==0 && ss(1)==0 && velRMS(step+1)<=sstol && accRMS(step+1)<=sstol && abs((thrustmag(step+1)- ...
            thrustmag(step))/thrustmag(step+1))/dt(step)<=sstol
        ss=[1,step];
    end
    if ss(1)==1
        break
    end

end
function [y]=simulateLimited2ndOrderActDiscrete(u_raw,tVec,w0, DampCoeff,rateLimit,lb,ub)

u_raw(u_raw>ub)=ub;
u_raw(u_raw<lb)=lb;
%setup
tStep=tVec(2)-tVec(1);


a = [0, 1; -w0^2, -2*DampCoeff*w0];
b = [0; w0^2];
c = [1, 0; 0, 1; -w0^2, -2*DampCoeff*w0];
d = [0; 0; w0^2];

dummyAct=c2d(ss(a,b,c,d),tStep);
a =dummyAct.a;
b =dummyAct.b;
c =dummyAct.c;
d =dummyAct.d;
%allocation
x=zeros(2,length(u_raw));
y=zeros(3,length(u_raw));

for iStep=1:length(tVec)-1
    
    x(:,iStep+1)=a*x(:,iStep)+b*u_raw(:,iStep);
    % rate limit
    if abs(x(2,iStep+1))>rateLimit
        x(2,iStep+1)=sign(x(2,iStep+1))*rateLimit;
        x(1,iStep+1)=(x(2,iStep+1)+x(2,iStep))/2*tStep+x(1,iStep);
    end
    % effective commanded def limit
    if (u_raw(:,iStep+1)-x(1,iStep+1))<-rateLimit*2*DampCoeff/w0
         u_raw(:,iStep+1)=-rateLimit*2*DampCoeff/w0+x(1,iStep+1);
    end
    % effective commanded def limit
    if (u_raw(:,iStep+1)-x(1,iStep+1))>rateLimit*2*DampCoeff/w0
        u_raw(:,iStep+1)=rateLimit*2*DampCoeff/w0+x(1,iStep+1);
    end
    % ub limit
    if x(1,iStep+1)>ub
        x(1,iStep+1)=ub;
        x(2,iStep)=(x(1,iStep+1)-x(1,iStep))/tStep;
    end
    % lb limit
    if x(1,iStep+1)<lb
        x(1,iStep+1)=lb;
        x(2,iStep)=(x(1,iStep+1)-x(1,iStep))/tStep;
    end
       
    y(:,iStep)=c*x(:,iStep)+d*u_raw(:,iStep);
end
%gradient repair
y(3,:)=gradient(gradient(y(1,:),tStep),tStep);
end
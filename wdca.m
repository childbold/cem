function [DH,DV] = wdca(R,phd,phpd,b0d,n,nu)

%Inputs
% r: distance from wedge
% phd: Diffraction angle (degrees)
% phpd: In1idence angle (degrees)
% b0d: beta_0 (degrees)
% n: Exterior wedge angle = n*pi
% nu: 1--> UTD, 2--> GTD
%   

sml=1e-8;
k0=2*pi;
ph=phd*pi/180;
php=phpd*pi/180;
sb0=sin(b0d*pi/180);

if nu==1			  %UTD
    if abs(ph-php+pi)<=sml %within face "n" ISB
        if(ph-php+pi)<0
            sgn = -1;
        else
            sgn =  1;
        end
		D1=n*sgn*sqrt(2*pi*k0*R)*exp(1i*pi/4);
		t2=((ph-php)-pi)/(2*n*pi);
        if t2>0
            nm1=floor(t2+0.5);
        else
            nm1=ceil(t2-0.5);
        end
		t3=(pi+(ph+php))/(2*n*pi);
        if t3>0
            np2=floor(t3+0.5);
        else
            np2=ceil(t3-0.5);
        end
        
        t4=((ph+php)-pi)/(2*n*pi);
        if t4>0
            nm2=floor(t4+0.5);
        else
            nm2=ceil(t4-0.5);
        end
		am1=1+cos(2*n*pi*nm1-(ph-php));
		ap2=1+cos(2*n*pi*np2-(ph+php));
		am2=1+cos(2*n*pi*nm2-(ph+php));
		X2=k0*R*am1;
		X3=k0*R*ap2;
		X4=k0*R*am2;
        [C2,S2]= fresnel(sqrt(2*X2/pi));
        [C3,S3]= fresnel(sqrt(2*X3/pi));
        [C4,S4]= fresnel(sqrt(2*X4/pi));
		F2=2*1i*sqrt(X2)*exp(1i*X2)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C2-1i*S2));
        F3=2*1i*sqrt(X3)*exp(1i*X3)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C3-1i*S3));        
        F4=2*1i*sqrt(X4)*exp(1i*X4)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C4-1i*S4));
		Cot2=cot((pi-(ph-php))/(2*n));
		Cot3=cot((pi+(ph+php))/(2*n));
		Cot4=cot((pi-(ph+php))/(2*n));
		D2=F2*Cot2;
		D3=F3*Cot3;
		D4=F4*Cot4;
    elseif abs(ph-php-pi)<=sml %within face "o" ISB       
        if(ph-php-pi) < 0
            sgn = -1;
        else
            sgn =  1;
        end
		D2=n*sgn*sqrt(2*pi*k0*R)*exp(1i*pi/4);
		t1=(pi+(ph-php))/(2*n*pi);        
        if t1>0
            np1=floor(t1+0.5);
        else
            np1=ceil(t1-0.5);
        end        
        t3=(pi+(ph+php))/(2*n*pi);        
        if t3>0
            np2=floor(t3+0.5);
        else
            np2=ceil(t3-0.5);
        end
		t4=((ph+php)-pi)/(2*n*pi);       
        if t4>0
            nm2=floor(t4+0.5);
        else
            nm2=ceil(t4-0.5);
        end
		ap1=1+cos(2*n*pi*np1-(ph-php));
		ap2=1+cos(2*n*pi*np2-(ph+php));
		am2=1+cos(2*n*pi*nm2-(ph+php));
		X1=k0*R*ap1;
		X3=k0*R*ap2;
		X4=k0*R*am2;
        [C1,S1]= fresnel(sqrt(2*X1/pi));
        [C3,S3]= fresnel(sqrt(2*X3/pi));
        [C4,S4]= fresnel(sqrt(2*X4/pi));
		F1=2*1i*sqrt(X1)*exp(1i*X1)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C1-1i*S1));
		F3=2*1i*sqrt(X3)*exp(1i*X3)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C3-1i*S3));
		F4=2*1i*sqrt(X4)*exp(1i*X4)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C4-1i*S4));
		Cot1=cot((pi+(ph-php))/(2*n));
		Cot3=cot((pi+(ph+php))/(2*n));
		Cot4=cot((pi-(ph+php))/(2*n));
		D1=F1*Cot1;
		D3=F3*Cot3;
		D4=F4*Cot4;
    elseif abs(ph+php-(2*n-1)*pi)<=sml  %within face "n" RSB
        if(ph+php-(2*n-1)*pi)<0
            sgn = -1;
        else
            sgn = 1;
        end
		D3=n*sgn*sqrt(2*pi*k0*R)*exp(1i*pi/4);
		t1=(pi+(ph-php))/(2*n*pi);
        if t1>0
            np1=floor(t1+0.5);
        else
            np1=ceil(t1-0.5);
        end
		t2=((ph-php)-pi)/(2*n*pi);
        if t2>0
            nm1=floor(t2+0.5);
        else
            nm1=ceil(t2-0.5);
        end
		t4=((ph+php)-pi)/(2*n*pi);
        if t4>0
            nm2=floor(t4+0.5);
        else
            nm2=ceil(t4-0.5);
        end
		ap1=1+cos(2*n*pi*np1-(ph-php));
		am1=1+cos(2*n*pi*nm1-(ph-php));
		am2=1+cos(2*n*pi*nm2-(ph+php));
		X1=k0*R*ap1;
		X2=k0*R*am1;
		X4=k0*R*am2;
        [C1,S1]= fresnel(sqrt(2*X1/pi));
        [C2,S2]= fresnel(sqrt(2*X2/pi));
        [C4,S4]= fresnel(sqrt(2*X4/pi));
		F1=2*1i*sqrt(X1)*exp(1i*X1)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C1-1i*S1));
		F2=2*1i*sqrt(X2)*exp(1i*X2)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C2-1i*S2));
		F4=2*1i*sqrt(X4)*exp(1i*X4)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C4-1i*S4));
		Cot1=cot((pi+(ph-php))/(2*n));
		Cot2=cot((pi-(ph-php))/(2*n));
		Cot4=cot((pi-(ph+php))/(2*n));
		D1=F1*Cot1;
		D2=F2*Cot2;
		D4=F4*Cot4;
    elseif abs(ph+php-pi)<=sml  %within face "o" RSB
        if (ph+php-pi) < 0.d0 
            sgn = -1;
        else
            sgn =  1;
        end
		D4=n*sgn*sqrt(2*pi*k0*R)*exp(1i*pi/4);
		t1=(pi+(ph-php))/(2*n*pi);
        if t1>0
            np1=floor(t1+0.5);
        else
            np1=ceil(t1-0.5);
        end
		t2=((ph-php)-pi)/(2*n*pi);
        if t2>0
            nm1=floor(t2+0.5);
        else
            nm1=ceil(t2-0.5);
        end
		t3=(pi+(ph+php))/(2*n*pi);
        if t3>0
            np2=floor(t3+0.5);
        else
            np2=ceil(t3-0.5);
        end
		ap1=1+cos(2*n*pi*np1-(ph-php));
		am1=1+cos(2*n*pi*nm1-(ph-php));
		ap2=1+cos(2*n*pi*np2-(ph+php));
		X1=k0*R*ap1;
		X2=k0*R*am1;
		X3=k0*R*ap2;
        [C1,S1]= fresnel(sqrt(2*X1/pi));
        [C2,S2]= fresnel(sqrt(2*X2/pi));
        [C3,S3]= fresnel(sqrt(2*X3/pi));
		F1=2*1i*sqrt(X1)*exp(1i*X1)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C1-1i*S1));
		F2=2*1i*sqrt(X2)*exp(1i*X2)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C2-1i*S2));
		F3=2*1i*sqrt(X3)*exp(1i*X3)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C3-1i*S3));
		Cot1=cot((pi+(ph-php))/(2*n));
		Cot2=cot((pi-(ph-php))/(2*n));
		Cot3=cot((pi+(ph+php))/(2*n));
		D1=F1*Cot1;
		D2=F2*Cot2;
		D3=F3*Cot3;
		else
		t1=(pi+(ph-php))/(2*n*pi);
        if t1>0
            np1=floor(t1+0.5);
        else
            np1=ceil(t1-0.5);
        end
		t2=((ph-php)-pi)/(2*n*pi);
        if t2>0
            nm1=floor(t2+0.5);
        else
            nm1=ceil(t2-0.5);
        end
		t3=(pi+(ph+php))/(2*n*pi);
        if t3>0
            np2=floor(t3+0.5);
        else
            np2=ceil(t3-0.5);
        end
		t4=((ph+php)-pi)/(2*n*pi);
        if t4>0
            nm2=floor(t4+0.5);
        else
            nm2=ceil(t4-0.5);
        end
		ap1=1+cos(2*n*pi*np1-(ph-php));
		am1=1+cos(2*n*pi*nm1-(ph-php));
		ap2=1+cos(2*n*pi*np2-(ph+php));
		am2=1+cos(2*n*pi*nm2-(ph+php));
		X1=k0*R*ap1;
		X2=k0*R*am1;
		X3=k0*R*ap2;
		X4=k0*R*am2;
        [C1,S1]= fresnel(sqrt(2*X1/pi));
        [C2,S2]= fresnel(sqrt(2*X2/pi));
        [C3,S3]= fresnel(sqrt(2*X3/pi));
        [C4,S4]= fresnel(sqrt(2*X4/pi));
		F1=2*1i*sqrt(X1)*exp(1i*X1)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C1-1i*S1));
		F2=2*1i*sqrt(X2)*exp(1i*X2)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C2-1i*S2));
		F3=2*1i*sqrt(X3)*exp(1i*X3)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C3-1i*S3));
		F4=2*1i*sqrt(X4)*exp(1i*X4)*sqrt(pi/2)*...
            ((exp(-1i*pi/4)/sqrt(2))-(C4-1i*S4));
		Cot1=cot((pi+(ph-php))/(2*n));
		Cot2=cot((pi-(ph-php))/(2*n));
		Cot3=cot((pi+(ph+php))/(2*n));
		Cot4=cot((pi-(ph+php))/(2*n));	
        D1=F1*Cot1;
		D2=F2*Cot2;
		D3=F3*Cot3;
		D4=F4*Cot4;
    end
    if (abs(php-pi)<sml)||(abs(php)<sml)   %grazing in1idence  
        DH=0;
        DV=-0.5*(exp(-1i*pi/4)/(2*n*sqrt(2*pi*k0)*sb0))...
                *((D1+D2)+(D3+D4));
    else
		DH=-(exp(-1i*pi/4)/(2*n*sqrt(2*pi*k0)*sb0))*((D1+D2)-(D3+D4));
		DV=-(exp(-1i*pi/4)/(2*n*sqrt(2*pi*k0)*sb0))*((D1+D2)+(D3+D4));
    end

	else		%GTD
		X=sin(pi/n)/(cos(pi/n)-cos((ph-php)/n));
		Y=sin(pi/n)/(cos(pi/n)-cos((ph+php)/n));
		DH=exp(-1i*pi/4)/(n*sqrt(2*pi*k0))*(X-Y);
		DV=exp(-1i*pi/4)/(n*sqrt(2*pi*k0))*(X+Y);
end
end




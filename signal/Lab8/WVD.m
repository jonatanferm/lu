function [WV,WV_smooth,A,A_smooth,k] = WVD(x,K,D,step,window,sigma)

% Calculate the analytical signal
    x = real(x);
    X = fft(x);
    if mod(length(x),2)
        Z = [1 2*ones(1,length(x)/2-1) zeros(1,length(x)/2)].*X;
    else
        Z = [1 2*ones(1,length(x)/2-1) 1 zeros(1,length(x)/2-1)].*X;
    end
    z = ifft(Z);
    clear Z X

% Define the window
    L=K+1;
    if window==1
        w=hamming(length(1:D:L));
        len=length(w);
        w=w(ceil(len/2):len);
    elseif window==0
        w=ones(length(1:D:K/2+1));
    end;

% Calculate the Wigner-Ville distribution
    z=[zeros(size(z,1),K/2) z zeros(size(z,1),K/2)];
    WV=[];
    j=0;
    for i=K/2+1:step:length(z)-K/2
        j=j+1;
        g1=w'.*z(i:D:i+K/2);
        g1=interp(g1,2);
        gk1=w'.*conj(z(i:-D:i-K/2));
        gk1=interp(gk1,2);
        h=(g1.*gk1)';
        hs=flipud(h);
        h=[hs(1:length(hs)-1);h];
        WV(:,j)=fft(h,2*K/D);
    end;
    if mod(j,2)
        WV = WV(:,1:end-1);
    end

% If wished, smooth the cross-terms
    if sigma
        % Transform to the ambiguity domain
            A = fft2(WV);
        % Define the kernel
            N = size(WV);
            nu = -N(1)/2:N(1)/2-1;
            tau = -N(2)/2:N(2)/2-1;
            k = exp(-nu(:).^2*tau.^2.*sigma);
        % Multiply in ambiguity domain and go back to Wigner-Ville domain
            A_smooth = A.*k;
            WV_smooth = ifft2(A_smooth);
            WV_smooth = WV_smooth(1:K/D,:);
    else
        k = [];
        A = [];
        A_smooth = [];
        WV_smooth = [];
    end
  
    WV = WV(1:K/D,:);

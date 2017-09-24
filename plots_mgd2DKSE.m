close all;
clear all;

fid = fopen('data_mgd2DKSE','r');
    output = textscan(fid,'%*s %f %*s %f %*s %f');
fclose(fid);

col1 = output{1}; col2 = output{2}; col3 = output{3};

M = col1(1); N = col2(1); dt = col3(1);
Tfinal = col1(2); p = col2(2); q = col3(2);
alpha = col1(3); beta = col2(3); gamma = col3(3);
sizeLvec = col1(4); L_1 = col2(4); stepL = col3(4);
sizedeltavec = col1(5); delta_1 = col2(5); stepdelta = col3(5);
a = col1(6);

L(1:sizeLvec)=0.0;
for i = 1:sizeLvec
    L(i) = L_1 + round(i-1)*stepL;
end

delta(1:sizedeltavec)=0.0;
for i = 1:sizedeltavec
    delta(i) = delta_1 + round(i-1)*stepdelta;
end

string_M = int2str(M);
string_N = int2str(N);

if M < 10
	string_M = ['000', string_M];
elseif M < 100
	string_M = ['00', string_M];
elseif M < 1000
	string_M = ['0', string_M];
elseif M < 10000
    string_M = string_M;
end

if N < 10
	string_N = ['000', string_N];
elseif N < 100
	string_N = ['00', string_N];
elseif N < 1000
	string_N = ['0', string_N];
elseif N < 10000
    string_N = string_N;
end

for j1 = 1:sizeLvec
    for j2 = 1:sizedeltavec
                  
Lint = floor(L(j1)); Lfrac = L(j1) - Lint;
aint = floor(a); afrac = a - aint;
deltaint = floor(delta(j2)); deltafrac = delta(j2) - deltaint;

Lint = num2str(Lint); aint = num2str(aint); deltaint = num2str(deltaint);

  Lfrac1 = floor(Lfrac*10);
  Lfrac2 = floor(Lfrac*100)-Lfrac1*10;
  Lfrac3 = floor(Lfrac*1000)-Lfrac1*100 - Lfrac2*10;
  
  afrac1 = floor(afrac*10);
  afrac2 = floor(afrac*100)-afrac1*10;
  afrac3 = floor(afrac*1000)-afrac1*100 - afrac2*10;
  
  deltafrac1 = floor(deltafrac*10);
  deltafrac2 = floor(deltafrac*100)-deltafrac1*10;
  deltafrac3 = floor(deltafrac*1000)-deltafrac1*100 - deltafrac2*10;

 Lfrac = [num2str(Lfrac1),num2str(Lfrac2),num2str(Lfrac3)];  

 afrac = [num2str(afrac1),num2str(afrac2),num2str(afrac3)];  

 deltafrac = [num2str(deltafrac1),num2str(deltafrac2),num2str(deltafrac3)];  

         % Load files:
         
         x_k1 = load(['x_k1_M=', string_M, '_N=', string_N, '.txt']);      
         y_k2 = load(['y_k2_M=', string_M, '_N=', string_N, '.txt']);  
         L2norminp = load(['L2norm_L=',Lint,'-',Lfrac,'_a=',aint,'-',afrac,'_delta=',deltaint,'-',deltafrac,'.txt']);        
         Linfnorminp = load(['Linfnorm_L=',Lint,'-',Lfrac,'_a=',aint,'-',afrac,'_delta=',deltaint,'-',deltafrac,'.txt']);
         profiles = load(['profiles_L=',Lint,'-',Lfrac,'_a=',aint,'-',afrac,'_delta=',deltaint,'-',deltafrac,'.txt']);     
         powerspecinp = load(['powerspec_L=',Lint,'-',Lfrac,'_a=',aint,'-',afrac,'_delta=',deltaint,'-',deltafrac,'.txt']);  
        
         x = x_k1(:,1); y = y_k2(:,1); k1 = x_k1(:,2); k2 = y_k2(:,2);
         uinitial = profiles(1:2*M,:); u = profiles(2*M+1:4*M,:); Fu_final = profiles(4*M+1:6*M,:);
         tdata = L2norminp(:,1); L2norm = L2norminp(:,2); Linfnorm = Linfnorminp(:,2);
         powerspec = powerspecinp(1:2*M,:);
 
         L1 = L(j1); L2 = L1^a ;
         
        % Plot of solution initial condition:        
        figure; 
        subplot(2,3,1)
        surfc(y*L2/L(end)^a,x*L1/L(end),uinitial);
        shading interp
        xlabel('$y$','Interpreter','LaTex','Fontsize',24); xlim([0 L2]);
        ylabel('$x$','Interpreter','LaTex','Fontsize',24); ylim([0 L1]);
        zlabel('$u_0$','Interpreter','LaTex','Fontsize',24);
        set(get(gca,'ZLabel'),'Rotation',0);
        
        % Plot of solution u:        
        %figure; 
        subplot(2,3,2)
        surfc(y*L2/L(end)^a,x*L1/L(end),u);
        shading interp
        xlabel('$y$','Interpreter','LaTex','Fontsize',24); xlim([0 L2]);
        ylabel('$x$','Interpreter','LaTex','Fontsize',24); ylim([0 L1]);
        zlabel('$u$','Interpreter','LaTex','Fontsize',24);
        set(get(gca,'ZLabel'),'Rotation',0);
        
        % Plot final spectrum Fu_final
        
        %figure; 
        subplot(2,3,3)
        surf(k2,k1,log10(Fu_final)); colorbar
        shading interp
          
        % Plot of L2norm squared:
        %figure; 
        subplot(2,3,4)
        plot(tdata,L2norm.^2,'-k');
        xlabel('t','Interpreter','LaTex','Fontsize',24);
        ylabel('$||u||_{2}^2$','Interpreter','LaTex','Fontsize',24);
        set(get(gca,'YLabel'),'Rotation',0);
      
        % Plot of Linfnorm:
        %figure; 
        subplot(2,3,5)
        plot(tdata,Linfnorm,'-k');
        xlabel('t','Interpreter','LaTex','Fontsize',24);
        ylabel('$||u||_{\infty}$','Interpreter','LaTex','Fontsize',24);
        set(get(gca,'YLabel'),'Rotation',0);
        
        % Plot powerspectrum
        %figure; 
        subplot(2,3,6)
        surf(k2,k1,log10(powerspec)); colorbar
        shading interp
     
		% Calculate mean of the L2norm squared over second half of time:
        nn = floor(max(size(L2norm))/2.0); mean=0;
        for i=1:nn
            mean = mean + L2norm(end-i+1)^2;
        end
        mean = mean/nn;
        fid = fopen('L2squaredmean.txt','a');
            fprintf(fid,'%f %f %f %f \n', L1, L2, a, mean);
        fclose(fid);
        
		% Calculate max of the Linfnorm over second half of time:
        nn = floor(max(size(L2norm))/2.0);
        maxval = max(Linfnorm(nn:end,1));
        fid = fopen('L2infmax.txt','a');
            fprintf(fid,'%f %f %f %f \n', L1, L2, a, maxval);
        fclose(fid);        
   
    end
end

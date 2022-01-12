tic
	clc;
	clear all;
	close all;
	load('C:\Users\Omistaja\Downloads\data\Part_1');
	FILE=[];
    
    %% mu own code for fininding the seven featurees:
    
    
	for d = 1:1
    Y = Part_1{1,d};
    p = Y(1, 1: 600);  %1min data
    bp = Y(2, 1:  600);
    e = Y(3,1: 600);
    figure(1)
    subplot(311)
    plot(p)
    xlabel('Number of Samples')
    ylabel('Amplitude')
    legend('PPG signal')
    title('Simple Overview of Collected PPG Signal')
    subplot(312)
    plot(bp)
    xlabel('Number of Samples')
    ylabel('Amplitude')
    legend('Truth BP signal')
    title('Simple Overview of Collected BP Signal')
    subplot(313)
    plot(e)
    legend('ECG signal')
    xlabel('Number of Samples')
    ylabel('Amplitude')
    title('Simple Overview of Collected ECG Signal')
    %% now find he second derivative of ppg
    p_1stD = gradient(p);
    figure (2)
    subplot(211)
    plot(p_1stD)
    p_2ndD = gradient(p_1stD);
    subplot(212)
    plot(p_2ndD)
    
    %% concatenation with 1 dim
    p_p2ndD_e = cat(1, p, p_2ndD, e);
    %plot(p_p2ndD_e)
    L=length(p_p2ndD_e);
	Ts=1/125; %sampling frequency=125Hz
	%t=linspace(0,10,0.008);
    T = 0: 1/125 : (L-1)/125 ;
    T1 = 0: 1/125 : (L-1)/125 ;
	%T =(0:0.008:7.999); %time vector based on sampling rate
    [pk,loc]= findpeaks(p); % max value of PPG signal
    figure(1)
    subplot(211)
    plot(T,p);
    hold on
    plot(T(loc),p(loc),'o')
    
    legend('PPG signal')
    xlabel('Time (Sec)')
    ylabel('Amplitude')
    title('Maximum Peak value detection for Systolic and Diastolic Time')
    
	PPG1=max(p)-p; % To find out the min peak of PPG
	[pk1,loc1]=findpeaks(PPG1,'MinPeakHeight',0.0); % min value of PPG signal
 
    subplot(212)
    plot(T1,PPG1);
    hold on
    plot(T1(loc1),PPG1(loc1),'o')
    legend('PPG signal')
    xlabel('Time (Sec)')
    ylabel('Amplitude')
    title('Minimum Peak value detection for Systolic and Diastolic Time')

%% find others parmeters

sys_time=0;
for i=1:1:length(loc)-1

    sys_time = sys_time + T1(loc(1,i))-T1(loc1(1,i));
end
sys_time = sys_time/length(loc)-1;



dias_time=0;
for i=1:1:length(loc)-1
    dias_time = dias_time + T(loc1(1,i))-T(loc(1,i));
end
plot(p)

dias_time = dias_time/length(loc)-1;

v = [0.1,0.25,0.33,0.5,0.66,0.75];

ppg_21_st = [];
ppg_21_dt = [];
for j=1:1:6
for i=loc1(1,1):1:loc(1,1)
    if(p(1,i)>=(v(j)*pk(1,1)+pk1(1,1)))
        a=i;
   
        break
    end
end

for i=loc(1,1):1:loc1(1,2)
  
    if(p(1,i)<=(v(j)*pk(1,1)+ pk1(1,1)))
        b=i;
        
        break
    end
end

ppg_21_st(j) = (loc(1,1)-a)*0.008;
ppg_21_dt(j) = (b-loc(1,1))*0.008;

end




%% 

    
    [pk2,loc2]= findpeaks(e,'MinPeakHeight',0.6); % max value of ECG signal
	%findpeaks(O1E,'MinPeakHeight',0.6);
	% original [pk3,loc3]=findpeaks(Fy1,'MinpeakHeight',0.0398); % max value of DPPG signal
    figure(3)
    subplot(211)
    plot(T1,e);
    hold on
    plot(T1(loc2),e(loc2),'o')
	[pk3,loc3]=findpeaks(p_2ndD,'MinPeakHeight',0.003); % max value of DPPG signal
    subplot(212)
    plot(T1,p_2ndD)
    hold on
    plot(T1(loc3), p_2ndD(loc3), 'o')
    
    
    %% finding the features
    [m,n] = size (loc2); % to find out vector dimensions of ECG signal
	[x,y] = size (loc3);
	
	P1=T(loc2); %ecg
	P=T(loc3); %derivative ppg
	P11=P1(1,1:n); %ecg
	P2= P(1,1:y); %derivative ppg
	ptt=0;
    temp=min(y,n);
	range=min(temp,5);
	for i=1:1:(range)
	
	ptt = ptt + abs(P2(1,(i))-P11(1,i));
	
	%PTT1(i) = P2(1,i)-P11(1,i) % To find out the transit time btwn ECG and PPG signal
    end
    
    ptt = ptt/(range);
    
    %% PIR calculation
    
    [lr,lr1] = size(loc1);  %min value of ppg
	rationum=0;
	ratioden=0;
	ih=0;
	il=0;
	for i=1:1:lr1-1
	rationum = rationum + pk(1,i);  % pk max value of PPG signal
	ratioden = ratioden + pk1(1,i);  % pk1 min value of PPG signal
	end
	%figure;
	
	ih = rationum/(lr1-1);
	il = ratioden/(lr1-1);
	
	PIR=ih/il;
    
    %% Hr calculation
    RR=diff(P1); % to find time taken for 1 heartbeat P1 = loc of ecg
	HR = 60./RR;
	hrfinal=0;
	[lr,lr1] = size (HR);
	tlr1 = lr1;
	for i=1:1:lr1
	t = HR(1,i);
	if t<=30||t>=200  %if condition full fill then less the number when make avg
	tlr1 = tlr1-1;
	else
	
	
	hrfinal = hrfinal + HR(1,i);
	end
	end
	hrfinal = hrfinal/(tlr1);

    
%% calculate the S, alpha and meu
% figure(1)
% 	subplot(3,1,1)
% 	plot(T,p);
% 	subplot(3,1,2)
% 	plot(T,e);
% 	subplot(3,1,3)
% 	plot(T,p_2ndD);
    
    
	Yy = fft(p);  %ppg signal
	% % P2 = abs(Y/L);
	%figure(2);%plot(Yy);
	Z=Yy(1);
	Yy(1)=0;
	S=real(ifft(Yy));
	
% 	figure(3);%plot(S);
	[pk4,loc4]=findpeaks(S); % max AC(ifft) value of PPG signal
% 	xlabel('Time(s)');
% 	ylabel('Ac amplitude(V)');
    
    [pk5,loc5]=findpeaks(bp);
	
	[lr,lr1] = size (loc4);
	iftmax=0;
	for i=1:1:lr1-1
	
	iftmax = iftmax + pk4(1,i);
	end
	
	meu = iftmax/(lr1-1);  %%  basically ppg signal , ifft peaks
	%figure;
	
	
	alpha = il*sqrt(1060*hrfinal/meu);   % il- ppg signal lowes value, hr-heart signal, 1060 for formula
    
    
    %% claculate the bp measurement
    %findpeaks(bp);
	BP1=max(bp)-bp; % To find out the min peak of BP
	[pk6,loc6]=findpeaks(BP1); % min value of BP(diastole) signal
	%findpeaks(BP1);
	
	
	[lr,lr1] = size (loc5);  %max bp signal loc
	bpmax=0;
	for i=1:1:lr1-1
	
	bpmax = bpmax + pk5(1,i);
	end
	
	bpmax = bpmax/(lr1-1);
	
	[lr,lr1] = size (loc6);  %min bp signal loc
	bpmin=0;
	for i=1:1:lr1-1
	
	bpmin = bpmin + pk6(1,i);
	end
	
	bpmin = bpmin/(lr1-1);
    
    
    
 %% all the parameters filed up
 filerow= [real(alpha) real(PIR) real(ptt) real(bpmax) real(bpmin) hrfinal ih il meu];
 filerow1= [ppg_21_dt(1) ppg_21_st(1)+ppg_21_dt(1) ppg_21_dt(1)/ppg_21_st(1) ppg_21_dt(2) ppg_21_st(2)+ppg_21_dt(2) ppg_21_dt(2)/ppg_21_st(2) ppg_21_dt(3) ppg_21_st(3)+ppg_21_dt(3) ppg_21_dt(3)/ppg_21_st(3) ppg_21_dt(4) ppg_21_st(4)+ppg_21_dt(4) ppg_21_dt(4)/ppg_21_st(4) ppg_21_dt(5) ppg_21_st(5)+ppg_21_dt(5) ppg_21_dt(5)/ppg_21_st(5) ppg_21_dt(6) ppg_21_st(6)+ppg_21_dt(6) ppg_21_dt(6)/ppg_21_st(6) sys_time dias_time]
 FILE = [FILE;filerow filerow1];
end
 
csvwrite('seven2polash.csv',FILE);
toc
 %iH-higest value of ppg peak, Il-lowest value of peak.  
 %PTT- pulse transit time between ECG and PPG.
 %PIR-Ratio between min and max value PPG.
 %meu = iftmax/(lr1-1);  %%  basically ppg signal , ifft peaks
 %alpha = il*sqrt(1060*hrfinal/meu);   % il- ppg signal lowes value, hr-heart signal, 1060 for formula

  
   
    
    
    
    
    
	%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tic
	clc;
	clear all;
	close all;
	load('C:\Users\Omistaja\Downloads\data\Part_1');
	FILE=[];
    
	for d=1:10
	d
	Y=(Part_1{1,d});
	O1P=Y(1,1:500);
	BP=Y(2,1:500);
	O1E=Y(3,1:500);
	%{
	figure('Name','ECG');
	plot(O1E);
	
	figure('Name','PPG');
	plot(O1P);
	%}
	[Fy]=gradient(O1P);
	%figure('Name','PPG 1st derivative');
	
	%plot(Fy);
	[Fy1]=gradient(Fy); % 2nd derivative
	%figure('Name', 'PPG 2nd derivative');
	
	%plot(Fy1);
	F=cat(1,O1P,Fy1,O1E);
	%figure('Name', 'All concat');
	%plot(F)
	L=length(F);
	Ts=1/125; %sampling frequency=125Hz
	%t=linspace(0,10,0.008);
	T =(0:0.008:7.999); %time vector based on sampling rate
	
	[pk,loc]= findpeaks(O1P); % max value of PPG signal
	PPG1=max(O1P)-O1P; % To find out the min peak of PPG
	[pk1,loc1]=findpeaks(PPG1,'MinPeakHeight',0.0); % min value of PPG signal
	
	%findpeaks(PPG1,'MinPeakHeight',0.6); % noise threshold
	%figure('Name','min PPG after threshold');
	
	%plot(PPG1);
	
	% original [pk2,loc2]= findpeaks(O1E,'MinpeakHeight',0.6); % max value of ECG signal
	
	[pk2,loc2]= findpeaks(O1E,'MinPeakHeight',0.6); % max value of ECG signal
	%findpeaks(O1E,'MinPeakHeight',0.6);
	% original [pk3,loc3]=findpeaks(Fy1,'MinpeakHeight',0.0398); % max value of DPPG signal
	[pk3,loc3]=findpeaks(Fy1,'MinPeakHeight',0.003); % max value of DPPG signal
	
	
	[m,n] = size (loc2); % to find out vector dimensions of ECG signal
	[x,y] = size (loc3);
	
	P1=T(loc2);
	P=T(loc3);
	P11=P1(1,1:n);
	P2= P(1,1:y);
	ptt=0;
	
	temp=min(y,n);
	range=min(temp,5);
	for i=1:1:(range)
	
	ptt = ptt + abs(P2(1,(i))-P11(1,i));
	
	%PTT1(i) = P2(1,i)-P11(1,i) % To find out the transit time btwn ECG and PPG signal
	end
	ptt = ptt/(range);
	
	[lr,lr1] = size(loc1);
	rationum=0;
	ratioden=0;
	ih=0;
	il=0;
	for i=1:1:lr1-1
	rationum = rationum + pk(1,i);
	ratioden = ratioden + pk1(1,i);
	end
	%figure;
	
	ih = rationum/(lr1-1);
	il = ratioden/(lr1-1);
	
	PIR=ih/il;
	
	RR=diff(P1); % to find time taken for 1 heartbeat
	HR = 60./RR;
	hrfinal=0;
	[lr,lr1] = size (HR);
	tlr1 = lr1;
	for i=1:1:lr1
	t = HR(1,i);
	if t<=30||t>=200
	tlr1 = tlr1-1;
	else
	
	
	hrfinal = hrfinal + HR(1,i);
	end
	end
	hrfinal = hrfinal/(tlr1);
	
	%figure
	%subplot(3,1,1)
	%plot(T,O1P);
	%subplot(3,1,2)
	%plot(T,O1E);
	%subplot(3,1,3)
	%plot(T,Fy1);
	Yy = fft(O1P);
	% % P2 = abs(Y/L);
	%figure;plot(Yy);
	Z=Yy(1);
	Yy(1)=0;
	S=real(ifft(Yy));
	
	%figure;plot(S);
	[pk4,loc4]=findpeaks(S); % max AC value of PPG signal
	xlabel('Time(s)');
	ylabel('Ac amplitude(V)');
	[pk5,loc5]=findpeaks(BP);
	
	[lr,lr1] = size (loc4);
	iftmax=0;
	for i=1:1:lr1-1
	
	iftmax = iftmax + pk4(1,i);
	end
	
	meu = iftmax/(lr1-1);
	%figure;
	
	
	alpha = il*sqrt(1060*hrfinal/meu);
	
	
	findpeaks(BP);
	BP1=max(BP)-BP; % To find out the min peak of BP
	[pk6,loc6]=findpeaks(BP1); % min value of BP(diastole) signal
	findpeaks(BP1);
	
	
	[lr,lr1] = size (loc5);
	bpmax=0;
	for i=1:1:lr1-1
	
	bpmax = bpmax + pk5(1,i);
	end
	
	bpmax = bpmax/(lr1-1);
	
	[lr,lr1] = size (loc6);
	bpmin=0;
	for i=1:1:lr1-1
	
	bpmin = bpmin + pk6(1,i);
	end
	
	bpmin = bpmin/(lr1-1);
	
	
	
	filerow= [real(alpha) real(PIR) real(ptt) real(bpmax) real(bpmin) hrfinal ih il meu];
	FILE = [FILE;filerow];
    end
    csvwrite('seven1polashdhar.csv',FILE);
	toc

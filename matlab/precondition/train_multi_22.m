 clear all
 close all
 format long
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %reading time domain data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fil = hdf5info('test.h5');
[a,b] = fil.GroupHierarchy.Groups(2).Datasets.Name;
c = fil.GroupHierarchy.Groups(11).Datasets(1).Name;
eob_A1{1} = h5read('test.h5' ,a );
eob_t{1} =h5read('test.h5', b);
eob_phi{1} = h5read('test.h5', c);

fil = hdf5info('test4096.h5');
[a,b] = fil.GroupHierarchy.Groups(2).Datasets.Name;
c = fil.GroupHierarchy.Groups(11).Datasets(1).Name;
eob_A1{2} = h5read('test4096.h5' ,a );
eob_t{2} =h5read('test4096.h5', b);
eob_phi{2} = h5read('test4096.h5', c);

fil = hdf5info('test3841.h5');
[a,b] = fil.GroupHierarchy.Groups(2).Datasets.Name;
c = fil.GroupHierarchy.Groups(11).Datasets(1).Name;
eob_A1{3} = h5read('test3841.h5' ,a );
eob_t{3} =h5read('test3841.h5', b);
eob_phi{3} = h5read('test3841.h5', c);

n=2; % choose how many waveform to be used max =3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Aligning waveform at t=0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:n; %M*M*M
    y = find(eob_A1{i} == max(eob_A1{i})); %gives index of max amplitude.

    eob_t1{i} = eob_t{i} - eob_t{i}(y); %shifting waveform t=0 corresponds to max(Amp)
    ini_time(i) = eob_t1{i}(1);
    end_time(i) =eob_t1{i}(end);
    %eob_phi{i} = unwrap(eob_phi{i});
end

end_t = min(end_time); 
ini_t = max(ini_time);

%choping off waveform at start and end according to smallest waveform
for i=1:n %M*M*M
    z = find(eob_t1{i} <= ini_t);
    z1 = find(eob_t1{i} >= end_t);
    eob_t1{i} = eob_t1{i}(z(end):z1(1));
    eob_A1{i} = eob_A1{i}(z(end):z1(1));
    eob_phi{i} = eob_phi{i}(z(end):z1(1));
end

for i =1:1:n; %M*M*M

    eob_h1{i} = eob_A1{i} .* exp(-1j .* eob_phi{i});
    eob_phi{i}= eob_phi{i} - eob_phi{i}(1); % setting initial phase zero
    t1 = find(eob_t1{i} <= -10^3); % indx where t=-10^3
    t(i) =t1(end);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %reduce data size, sampling in phase
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:n; %M*M*M

    %phas = wrapToPi(eob_phi{k});
    
   
     t1   = eob_t1{k}(1:t(k)-1);
     A1    = eob_A1{k}(1:t(k)-1);
     p1    = (eob_phi{k}(1:t(k)-1));
     
% resample: uniform sampling in freq (interp)
% %  N = 5000;
% % f1 = diff(p1)./diff(t1);
% % f2 = [f1(1): (f1(end)-f1(1)) /N:f1(end)];
% t3 = t1(1):100:t1(end);
% time1 =t3;
% time1 = spline(t1, t1(1:end), t3);
% phase1  = spline(t1, p1(1:end), t3);
% amp1 = spline(t1, A1(1:end), t3);

% resample: in phase also keeping peak and lowest point from wraptoPi phase
p11 = diff(wrapToPi(p1));
pindx = find((p11) <= -6.1);

l = [pindx];
for j=1:length(l)
    l33(j) = l(j)+1;
end
m =1;
for j = p1(1):pi:p1(end)
    l11 = find(p1 <= j);
    l1(m) = l11(end);
    m =m+1;
    clear l11;
end
l2 = [l1,l',l33];

[b,m1,n1] = unique(l2,'first');
[c1,d1] =sort(m1);
b = b(d1);

l3 = sort(b);
pu = (p1(l3))';


% pu = [pu, p1(end)];
% %pu = [p1(1):pi/4:p1(end)];
%pu = p1(1):pi:p1(end);

phase1 = (pu);
time1 = spline(p1, t1, pu);
%phase1  = spline(p1, p1, pu);
amp1 = spline(p1, A1, pu);
    
    % uniform in time
    
     t2   = eob_t1{k}(t(k)-1:end);
     A2    = eob_A1{k}(t(k)-1:end);
     p2    = eob_phi{k}(t(k)-1:end);
     
     dt=0.1;
     time2 = t2(1)+0.1:dt:t2(end);
     
    phase2 = spline(t2, p2, time2);
    amp2 = spline(t2, A2, time2);
    
    % combining both
    eob_AA{k} = [amp1, amp2];
    eob_phii{k} = [phase1, phase2];
    eob_time{k} = [time1, time2];
    tt(k) = find(eob_time{k} == time1(end));
    
   clear amp1; clear amp2; clear phase1; clear phase2; clear time1; clear time2; clear A1; clear A2; clear t1; clear t2; clear p1;
   clear p2; clear l; clear l1; clear pu; clear l3; clear p11; clear pindx; 
    
 end
 %Now choose the smallest waveform
len = cellfun(@length, eob_time);
small_wf = find( len == min(len));

% constructing h

for i =1:1:n; %M*M*M
    h{i} = eob_AA{i} .* exp(-1j .* eob_phii{i});

    
end

% interpolating on smallest grid
 for i =1:1:n; %M*M*M
     eob_pha4{i}= spline( eob_time{i}, eob_phii{i}, eob_time{small_wf});
     
     eob_A4{i}= spline( eob_time{i}, eob_AA{i}, eob_time{small_wf});
    
 end
    


% constructing h

for i =1:1:n; %M*M*M

    eob_h{i} = eob_A4{i} .* exp(-1j .* eob_pha4{i});
end




inst_p=angle((h{1}));
inst_p1=angle((eob_h{1}));

% plotting

% plot(eob_time{1}, inst_p1)
% hold on
% plot(eob_time{1}, inst_p)
% legend('before','after interpolated')
% title('instantaneous phase')

figure(1)
plot(eob_time{small_wf}, eob_h{1})
title('interpolated waveform')
% 
%  figure(3)
% subplot(3,1,2)
% plot(eob_time{1},eob_AA{1} )
% title('Amplitude -data resampled & shifted')
% 
% subplot(3,1, 3)
% plot(eob_time{2},eob_A{1} )
% title('interpolated and shifted')
% 
%  
%  
%  figure(2)
%  subplot(3,1,2)
% plot(eob_time{1},eob_phii{1} )
% title('phase -data resampled & shifted')
% 
% subplot(3,1, 3)
% plot(eob_time{2},eob_pha{1} )
% title('interpolated and shifted')
%  subplot(3,1,1)
%  plot(eob_t1{1},eob_phi{1})
%  title('data-shifted')
 

 


%%

clear;clc
cd('E:\LeeLab-Coop1\Frozen tissue\sample12')
load('circle_ave_intensity_S_new')
p1=[535:540 565:585];%AD1769
p2=[28 29 30 31];%AD1769
%cd('G:\23IDB_2015_08_07\isolate\raster2')
%load('circle_ave_intensity_S_new')
%p1=[117:123 158:164 199:205];%isolate
%p2=[1:9];%isolate
%cd('G:\23IDB_2015_08_07\picks2\raster3')
%load('circle_ave_intensity_S_new')
%p1=[1:5 11:15 21:25];%picks
%p2=[87:90 97:100];%picks
%cd('G:\23IDB_2015_08_07\mismatch2\raster1')
%load('circle_ave_intensity_S_new')
%p1=[2401:2405 2451:2455];%mismatch2_rast1
%p2=[1397:1400 1447:1450];%mismatch2_rast1
%cd('G:\23IDB_2015_08_07\mismatch2\raster2')
%load('circle_ave_intensity_S_new')
%p1=[81:85 91:95];%mismatch2_rast2
%p2=[8:10 18:20];%mismatch2_rast2
%cd('G:\23IDB_2015_08_07\AD3\raster1')
%load('circle_ave_intensity_S_new')
%p1=[824:828 864:868];%AD3
%p2=[540:545];%AD3

IB=sum(I(:,p1)')'/length(p1);
IB1=smooth(sum(I(:,p2)'/length(p2))'-IB,0.01,'rloess');%intensity of tissues, used as control%%

A=[IB(100:720)./sum(IB(100:720)) IB1(100:720)./sum(IB1(100:720))]';%% august sample setting
cd('E:\LeeLab-Coop1\Frozen tissue\sample12')
load('circle_ave_intensity_S_new')
for i=1:length(I(1,:))
    %figure(2);hold on;plot(r1(30:end),I(30:end,i)./sum(I(30:end,i)))
    b=I(100:720,i)./sum(I(100:720,i));
    cvx_begin
        variable x(2)
        minimize(norm(b'-x'*A))
        subject to
             x(1) >= 0
             x(2) >= 0
             x(1)+x(2) == 1

    cvx_end
    para(:,i) = x;
    I_tissue(30:900,i)=I(30:900,i)./sum(I(30:900,i))-x(1)*IB(30:900)./sum(IB(30:900));
end
save('new_sub_bkgd1','I_tissue','r1','I','IB','IB1','para')
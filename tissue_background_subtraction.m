clear,clc
cd('E:\LeeLab-Coop1\test 6 frozen')
load('new_sub_bkgd1')
New1 = [24];
for m = 1
    New = New1(m,:);
    I2=zeros(1,900);
    I2(30:900)=(sum(I_tissue(30:900,New)')/length(New));
end

cd('E:\LeeLab-Coop1\test 6 frozen')
load('new_sub_bkgd1')
clear('para')
A = [I2; ones(1,900)];
%sample 6 frozen parameters
A1 = [A(:,400:520) A(:,550:570)];
%A1 = [A(:,300:700)];
A2 = A(:,30:900);
for i = 1:960
    I_tissue(30:900,i) = smooth(I_tissue(30:900,i),11);
    b = I_tissue(:,i);
    b = smooth(b,20,'rloess');
    %sample 6 frozen parameters
    b = [b(370:490);b(520:540)];
    %b = [b(300:700)];
    cvx_begin
        variable x(2)
        minimize(norm(b'-x'*A1))
        subject to
                min(x) >= 0
    cvx_end
    para(1:2,i) = x;
    I3(30:900,i)=(I_tissue(30:900,i)'-x'*A2).';
end
plot(I3)
save('further_bkgd_procs','I3','para','A2','I_tissue','r1')

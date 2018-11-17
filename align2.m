

function [stk2]=align2(stk,showim,iteration) %stk is an image stack
 
[y x z] = size(stk); 
stk(stk<0)=0;
n=2; step=1; % n = number of steps, step = step size (1pix)
stk2=zeros(y,x,z); stk2(:,:,1)=stk(:,:,1);
sortarray1=zeros(2,z); sortarray2=zeros(2,z);
corr=zeros(2*n+1,2*n+1,z);
 
%% Normalization
for mm=1:z
    stk(:,:,mm)=stk(:,:,mm)/max(max(stk(:,:,mm)));
end
reference=stk(:,:,1);  %reference (1st) image
 
%% Cross-Correlation normal orientation
for mm=2:z
    for jj=-n:n
        for kk=-n:n
            shift=squeeze(circshift(stk(:,:,mm),[jj*step kk*step]));
            shift(shift<0.2)=0;
            diff=abs(reference(n*step+1:y-n*step,n*step+1:x-n*step)-shift(n*step+1:y-n*step,n*step+1:x-n*step));
            corr(jj+(n+1),kk+(n+1),mm)=(sum(sum(diff)));
        end
    end
    MIN1=min(min(corr(:,:,mm)));
    [S1,S2]=find(corr(:,:,mm)==MIN1);
    sortarray1(1,mm)=(S1(1)-(n+1))*step; sortarray1(2,mm)=(S2(1)-(n+1))*step;
end
 
%% Cross-Correlation flipped 180 degrees
stktemp=stk;
for mm=2:z
    if mm>1
        stktemp(:,:,mm)=rot90(stktemp(:,:,mm),2); 
    end
    for jj=-n:n
        for kk=-n:n
            shift=squeeze(circshift(stktemp(:,:,mm),[jj*step kk*step]));
            shift(shift<0.2)=0;
            diff=abs(reference(n*step+1:y-n*step,n*step+1:x-n*step)-shift(n*step+1:y-n*step,n*step+1:x-n*step));
            corr(jj+(n+1),kk+(n+1),mm)=(sum(sum(diff)));
        end
    end
    MIN2=min(min(corr(:,:,mm)));
    [S1,S2]=find(corr(:,:,mm)==MIN2);
    sortarray2(1,mm)=(S1(1)-(n+1))*step; sortarray2(2,mm)=(S2(1)-(n+1))*step;
end
 
%% Realing images
if MIN1<MIN2
    for mm=2:z
        S1=sortarray1(1,mm); S2=sortarray1(2,mm);
        stk2(:,:,mm)=squeeze(circshift(stk(:,:,mm),[S1 S2])); 
    end
else
    for mm=2:z
        S1=sortarray2(1,mm); S2=sortarray2(2,mm);
        stk2(:,:,mm)=squeeze(circshift(stktemp(:,:,mm),[S1 S2])); 
    end
end
 
%% Figures
if showim==1
figure(101), 
subplot(1,3,1), imagesc(squeeze(stk(:,:,1))), axis image, title(strcat('reference img',int2str(iteration)));
subplot(1,3,2), imagesc(squeeze(stk(:,:,2))), axis image, title 'new img';
subplot(1,3,3), imagesc(squeeze(stk2(:,:,2))), axis image, title 'shifted img';
end

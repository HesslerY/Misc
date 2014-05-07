load sol.mat

plot(Ap(:,3))

for i=1:length(Ap(:,1))-2
    if abs(Ap(i+1,3)-Ap(i,3))/abs(Ap(i,3))>0.05
    Ap(i+1,3)=NaN;
    end
    if abs(Ap(i+2,3)-Ap(i,3))/abs(Ap(i,3))>0.05
    Ap(i+2,3)=NaN;
    end
end

nanmean(Ap(:,3))


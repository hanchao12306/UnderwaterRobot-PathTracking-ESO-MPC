function step = step_cal(Pobs,len_step,D,Pi)
%NEAROBS 此处显示有关此函数的摘要
%   此处显示详细说明

step_max=4*len_step;
step_min=0.8*len_step;
sum=0;
nummax=6;
nummin=1;
for k=1:length(Pobs(:,1))
    dist=(Pobs(k,1)-Pi(1))^2+(Pobs(k,2)-Pi(2))^2;
    if dist<D^2
       sum=sum+1; 
    end
end
if sum>nummax || sum==nummax
   step=step_min; 
end
if sum<nummin || sum==nummin
   step=step_max; 
end
if sum<nummax && sum>nummin
    step=step_max-(sum-nummin)*(step_max-step_min)/(nummax-nummin);
end


end


function res_c=get_residual_C(c)
% input: Cb or Cr
res_c=zeros(size(c));
pred_c=zeros(size(c));
[m,n]=size(c);
res_c(1,:)=c(1,:);
res_c(:,1)=c(:,1);
pred_c(1,:)=c(1,:);
pred_c(:,1)=c(:,1);
for i=2:m
    for j=2:n
        pred_c(i,j)=(3/8)*pred_c(i,j-1)-(1/4)*pred_c(i-1,j-1)+(7/8)*pred_c(i-1,j);
        res_c(i,j)=round(c(i,j)-pred_c(i,j));
        pred_c(i,j)=pred_c(i,j)+res_c(i,j);
    end
end   
end
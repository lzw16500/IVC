function res_y=get_residual_Y(y)
%input: Y
res_y=zeros(size(y));
pred_y=zeros(size(y));
[m,n]=size(y);
res_y(1,:)=y(1,:);
res_y(:,1)=y(:,1);
pred_y(1,:)=y(1,:);
pred_y(:,1)=y(:,1);
for i=2:m
    for j=2:n
        pred_y(i,j)=(7/8)*pred_y(i,j-1)-(1/2)*pred_y(i-1,j-1)+(5/8)*pred_y(i-1,j);
        res_y(i,j)=round(y(i,j)-pred_y(i,j));
        pred_y(i,j)=pred_y(i,j)+res_y(i,j);
    end
end    
end
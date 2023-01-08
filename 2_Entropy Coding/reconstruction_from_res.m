function [rec_y,rec_down_cb,rec_down_cr]=reconstruction_from_res(res_y,res_cb,res_cr) % reconstruct the image based on the residual matrix
[m1,n1]=size(res_y);
[m2,n2]=size(res_cb);
rec_y=zeros(m1,n1);
rec_down_cb=zeros(m2,n2);
rec_down_cr=zeros(m2,n2);
rec_y(1,:)=res_y(1,:);
rec_down_cb(1,:)=res_cb(1,:);
rec_down_cr(1,:)=res_cr(1,:);
rec_y(:,1)=res_y(:,1);
rec_down_cb(:,1)=res_cb(:,1);
rec_down_cr(:,1)=res_cr(:,1);

for i=2:m1
    for j=2:n1
        rec_y(i,j)=(7/8)*rec_y(i,j-1)-(1/2)*rec_y(i-1,j-1)+(5/8)*rec_y(i-1,j)+res_y(i,j);
    end
end

for i=2:m2
    for j=2:n2
        rec_down_cb(i,j)=(3/8)*rec_down_cb(i,j-1)-(1/4)*rec_down_cb(i-1,j-1)+(7/8)*rec_down_cb(i-1,j)+res_cb(i,j);
        rec_down_cr(i,j)=(3/8)*rec_down_cr(i,j-1)-(1/4)*rec_down_cr(i-1,j-1)+(7/8)*rec_down_cr(i-1,j)+res_cr(i,j);
    end
end

end
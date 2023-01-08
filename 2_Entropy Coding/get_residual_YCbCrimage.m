function res_im= get_residual_YCbCrimage(image) % res_im has the same size as image
[m,n,c]=size(image);
res_im=zeros(m,n,c);
predict_im=zeros(m,n,c);
res_im(:,1,:)=image(:,1,:);
res_im(1,:,:)=image(1,:,:);
predict_im(:,1,:)=image(:,1,:);
predict_im(1,:,:)=image(1,:,:);
for i=1:c
    if i==1
        for j=2:m
            for k=2:n
                predict_im(j,k,i)=(7/8)*predict_im(j,k-1,i)-(1/2)*predict_im(j-1,k-1,i)+(5/8)*predict_im(j-1,k,i);
                res_im(j,k,i)=round(image(j,k,i)-predict_im(j,k,i));
                predict_im(j,k,i)=predict_im(j,k,i)+res_im(j,k,i);
            end
        end
    else
        for j=2:m
            for k=2:n
                predict_im(j,k,i)=(3/8)*predict_im(j,k-1,i)-(1/4)*predict_im(j-1,k-1,i)+(7/8)*predict_im(j-1,k,i);
                res_im(j,k,i)=round(image(j,k,i)-predict_im(j,k,i));
                predict_im(j,k,i)=predict_im(j,k,i)+res_im(j,k,i);
            end
        end
    end
end

end
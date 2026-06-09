clc;
clear;

load('region_data.mat');

data_dims = size(region_data);
num_times = data_dims(1);
num_features = data_dims(2);

% CR
cr_data = region_data(45:num_times, 1, :, :);  
data_dims_cr = size(cr_data);
num_times_cr = data_dims_cr(1);
cr_matrix = reshape(cr_data, num_times_cr, 16*16);  

mydata=cr_matrix;
win_m=10;

lamb=0;
    lamb=lamb+1;
    flag=1;
    for ii=1:size(mydata,1)-win_m+1 
        xx=mydata(ii:size(mydata,1),:)';
        noisestrength=0;          % noise could be added

        xx_noise=xx+noisestrength*rand(size(xx));
        [input_dimensions, time_points]=size(xx);
        trainlength=win_m;

        flat_y=zeros(trainlength,1);

        embeddings_num=4;   

        traindata=xx_noise(:,1:trainlength);
        traindata=traindata-mean(traindata,2);
      
        X=traindata';       
        m=trainlength;
        n=input_dimensions;
        L=embeddings_num;
        P=X(2:end,:);
        Q=X(1:end-1,:);

      
        a=0.1;
        b=1-a;

        %%  solve Z  %%
        H=zeros(n*L,n*L);
        H(1:n, 1:n)=a*X'*X-b*P'*P;
        H(1:n, n+1:2*n)=b*P'*Q;

        for j=2:L-1
            H(n*(j-1)+1:n*j, n*(j-1)+1-n:n*j-n)=b*Q'*P;
            H(n*(j-1)+1:n*j, n*(j-1)+1:n*j)=a*X'*X-b*P'*P-b*Q'*Q;
            H(n*(j-1)+1:n*j, n*(j-1)+1+n:n*j+n)=b*P'*Q;
        end

        H(n*(L-1)+1:n*L, n*(L-2)+1:n*(L-1))=b*Q'*P;
        H(n*(L-1)+1:n*L, n*(L-1)+1:n*L)=a*X'*X-b*Q'*Q;

        [V,D]=eig(H);

 [weight_W_1,xi]=fpca_rebuild(H,n,L);


Z1=X*weight_W_1;
        if ii==1
            weight_Z=norm(Z1,'fro');
        end

        %  flat Z
        clear flat_z_1 flat_z_pred_1
        sd_flat_z_1=0;
        sd_flat_z_pred_1=0;
        for z1i=1:size(Z1,1)
            num=0;
            for z1j=1:size(Z1,2)
                if z1i-z1j+1<1
                    break;
                end
                num=num+1;
                tmp(num)=Z1(z1i-z1j+1,z1j);
              
            end
            tmp=tmp(1:num);
            flat_z_1(z1i)=mean(tmp);
            sd_flat_z_1=sd_flat_z_1+std(tmp);
        end
        sd_flat_z_1=sd_flat_z_1/size(Z1,1);

        for z1j=2:size(Z1,2)
            num=size(Z1,2)-z1j+1;
            for ni=1:num
                z1i=size(Z1,1)-ni+1;
                tmp(ni)=Z1(z1i,z1j+ni-1);
            end
            tmp=tmp(1:num);
            flat_z_pred_1(z1j-1)=mean(tmp);
            sd_flat_z_pred_1=sd_flat_z_pred_1+std(tmp);
            %             num
        end
        sd_flat_z_pred_1=sd_flat_z_pred_1/(size(Z1,2)-1);


        temp_var_y_1(ii)=norm(Z1,'fro');
        all_var_y_1(ii,lamb)= temp_var_y_1(ii);

        temp_var_flat_y_1(ii)=std(flat_z_1);
        temp_var_flat_y_pred_1(ii)=std(flat_z_pred_1);


        all_flat_z_1(ii,:)=flat_z_1;
        all_flat_z_pred_1(ii,:)=flat_z_pred_1;

    end

figure;
plot(temp_var_y_1, 'g*-','LineWidth',2);
title('SD of Z1');
set(gca,'FontSize',20);



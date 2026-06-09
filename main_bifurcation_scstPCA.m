clc;
warning off;
clear;
close all;

load sol.mat;
mydata=sol(:,2000:10:20000)';
win_m=8;


lamb=0;
    lamb=lamb+1;
    flag=1;
    clear temp_var_y
    for ii=1:size(mydata,1)-win_m+1 
        xx=mydata(ii:size(mydata,1),:)';
        noisestrength=0;          % noise could be added

        xx_noise=xx+noisestrength*rand(size(xx));
        [input_dimensions, time_points]=size(xx);
        trainlength=win_m;
        flat_y=zeros(trainlength,1);

        embeddings_num=3;    

        traindata=xx_noise(:,1:trainlength);
        traindata=traindata-mean(traindata,2);
   
        X=traindata';    

        m=trainlength;
        n=input_dimensions;
        L=embeddings_num;
        P=X(2:end,:);
        Q=X(1:end-1,:);

      
        a=0.001;
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

        ao=(diag(D));
        ao=real(ao);
        [aa, eigvIdx]=sort(ao,'descend');
        V=V(:,eigvIdx);

 [weight_W_1,xi]=fpca_rebuild(H,n,L);
all_xi(1:length(xi),ii,lamb)=xi/sum(xi);
weight_W_1=weight_W_1/norm(weight_W_1,'fro');
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
                %                 [num2str(zi-zj+1),', ', num2str(zj)]
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
                %                  [num2str(zi),', ', num2str(zj+ni-1)]
            end
            tmp=tmp(1:num);
            flat_z_pred_1(z1j-1)=mean(tmp);
            sd_flat_z_pred_1=sd_flat_z_pred_1+std(tmp);
            %             num
        end
        sd_flat_z_pred_1=sd_flat_z_pred_1/(size(Z1,2)-1);


        temp_var_y_1(ii)=norm(Z1,'fro');
        all_var_y_1(ii,lamb)= temp_var_y_1(ii);
        % all_embedding_error_1(ii,lamb)=embedding_error(Z1);

        temp_var_flat_y_1(ii)=std(flat_z_1);
        temp_var_flat_y_pred_1(ii)=std(flat_z_pred_1);


        all_flat_z_1(ii,:)=flat_z_1;
        all_flat_z_pred_1(ii,:)=flat_z_pred_1;

    end

%% BOCD(Rbest）
out = beast(temp_var_flat_y_1, ...
            'season', 'none', ...
            'tcp.minmax', [0, 1], ...          
            'torder.minmax', [0, 1], ...       
            'tseg.min', 10, ...                 
            'tseg.leftmargin', 10, ...         
            'tseg.rightmargin', 10);            
disp(out.trend.ncpPr);  

plotbeast(out);
prob_cp = out.trend.cpOccPr;   
prob_loc = out.trend.cp;


[max_prob,loc] = max(prob_cp);
fprintf('Most likely change point position: point %d, probability = %.4f\n', loc, max_prob);

figure;
plot(temp_var_y_1, '*-', 'LineWidth', 2, 'Color',[147,110,190] / 255);%126, 47, 142
title('SD of Z1');
set(gca,'FontSize',20);
hold on;  
yl = ylim;  
line([loc, loc], yl+[0,0.2], 'Color', [255, 51, 51]/255, 'LineStyle', '--', 'LineWidth', 3);
hold off;
clc;
warning off;
clear;
close all;

load filled_data.mat;
filled_data=filled_data(200:289,:);%155:219   213:289   200:289
data=filled_data(:, [3,4,5,6,7,8]);

mydata= flipud(data);
win_m=8;

lamb=0;
    lamb=lamb+1;
    flag=1;
    clear temp_var_y
    for ii=1:size(mydata,1)-win_m+1 
        xx=mydata(ii:size(mydata,1),:)';
        noisestrength=0;         

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

        ao=(diag(D));
        ao=real(ao);
        [aa, eigvIdx]=sort(ao,'descend');
        V=V(:,eigvIdx);

        for ci=1:length(aa)
            if aa(ci)>0
                break;
            end
        end
        eigv(ii)=aa(ci);
        cW=V(:,ci);

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
        end
        sd_flat_z_pred_1=sd_flat_z_pred_1/(size(Z1,2)-1);


        temp_var_y_1(ii)=norm(Z1,'fro');
        all_var_y_1(ii,lamb)= temp_var_y_1(ii);

        temp_var_flat_y_1(ii)=std(flat_z_1);
        temp_var_flat_y_pred_1(ii)=std(flat_z_pred_1);


        all_flat_z_1(ii,:)=flat_z_1;
        all_flat_z_pred_1(ii,:)=flat_z_pred_1;

    end
% end
age=flipud(filled_data(end-length(temp_var_y)-1:end,2));
out = beast(temp_var_flat_y_1, ...
            'season', 'none', ...
            'tcp.minmax', [0, 5], ...         
            'torder.minmax', [0, 2], ...       
            'tseg.min', 10, ...                  
            'tseg.leftmargin', 10, ...          
            'tseg.rightmargin', 10);             

cp_positions = out.trend.cp;    
cp_probs = out.trend.cpPr;        

valid = ~isnan(out.trend.cp);
cp_pos = out.trend.cp(valid);
cp_prob = out.trend.cpPr(valid);

thresh = 0.5;
keep = cp_prob > thresh;
cp_final = cp_pos(keep);
prob_final = cp_prob(keep);

cp_time = age(cp_final);

for i = 1:length(cp_final)
    fprintf('Change point %d：index %d', i, cp_final(i));
    fprintf(' (time %.2f)', cp_time(i));
    fprintf('，probability %.4f\n', prob_final(i));
end

figure('Position', [100, 100, 1200, 800]);

subplot(3,1,1);
plot(filled_data(:,2), filled_data(:,3), 'g-', 'LineWidth', 1);
xlabel('Age (cal yr BP)');
ylabel('calcite (%)');
title('Percentage of calcite');
grid on;
set(gca, 'XDir', 'reverse'); 

subplot(3,1,2);
plot(filled_data(:,2), filled_data(:,4), 'r-', 'LineWidth', 1);
xlabel('Age (cal yr BP)');
ylabel('δ¹³C (‰ VPDB)');
title('Carbon Isotope Record');
grid on;
set(gca, 'XDir', 'reverse');

subplot(3,1,3);
plot(filled_data(:,2), filled_data(:,5), 'b-', 'LineWidth', 1);
xlabel('Age (cal yr BP)');
ylabel('δ¹⁸O (‰ VPDB)');
title('Oxygen Isotope Record');
grid on;
set(gca, 'XDir', 'reverse'); 

figure('Position', [100, 100, 1200, 800]);

subplot(3,1,1);
plot(filled_data(:,2), filled_data(:,6), 'r-', 'LineWidth', 1);
xlabel('Age (cal yr BP)');
ylabel('δ¹³C_ac (‰ VPDB)');
title('Corrected Carbon Isotope Record');
grid on;
set(gca, 'XDir', 'reverse');

subplot(3,1,2);
plot(filled_data(:,2), filled_data(:,7), 'b-', 'LineWidth', 1);
xlabel('Age (cal yr BP)');
ylabel('δ¹⁸O_ac (‰ VPDB)');
title('Corrected Oxygen Isotope Record');
grid on;
set(gca, 'XDir', 'reverse');

subplot(3,1,3);
%valid_growth = ~isnan(data_clean.growth_mmyr);
plot(filled_data(:,2), filled_data(:,8), 'c-', 'LineWidth', 1);
xlabel('Age (cal yr BP)');
ylabel('Growth Rate (mm/yr)');
title('Stalagmite Growth Rate');
grid on;
set(gca, 'XDir', 'reverse');

age=flipud(filled_data(:,2))';
figure;
plot(age(1:length(temp_var_y)),temp_var_y, '*-','LineWidth',2);

title('SD of Z');
set(gca,'FontSize',20);
set(gca, 'XDir', 'reverse');

figure;
plot(age(1:length(temp_var_y)),temp_var_y_1, '*-','LineWidth',2,'Color',[112,206,188] / 255);%组合W
hold on;
xline(cp_time, 'r--', 'LineWidth', 2); 
hold off;
title('SD of Z1');
set(gca,'FontSize',20);
set(gca, 'XDir', 'reverse');


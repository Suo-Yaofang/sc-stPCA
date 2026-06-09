clc;
clear;
close all;
Y=mylorenz(10,0.02);
noisestrength=0;
XX=Y+noisestrength*rand(size(Y));% noise could be added

start=1;
start=randi([1 1000]);
ii=start-1;
num_rebuild_nice1=0;
for ii=155:1154
     ii=ii+1;
    xx=XX(3000+ii:size(XX,1), :)';
    [input_dimensions, time_point]=size(xx);
    
    noisestrength=0;          % noise could be added

    xx_noise=xx+noisestrength*rand(size(xx));
    trainlength=8;
    flat_y=zeros(trainlength,1);
    
    embeddings_num=3;    % dimensions after dimension reduction
    
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
    ci
    aa(ci)

    cW=V(:,ci);
    W=reshape(cW, n, L);
    Z=X*W;

    weight_W=W;
 
    if ii==start
        weight_Z=norm(Z,'fro');
    end

    ii

    clear flat_z flat_z_pred
    sd_flat_z=0;
    sd_flat_z_pred=0;
    for zi=1:size(Z,1)
        num=0;
        for zj=1:size(Z,2)
            if zi-zj+1<1
                break;
            end
            num=num+1;
            tmp(num)=Z(zi-zj+1,zj);
        end
        tmp=tmp(1:num);
        flat_z(zi)=mean(tmp);
        sd_flat_z=sd_flat_z+std(tmp);
    end
    sd_flat_z=sd_flat_z/size(Z,1);

    for zj=2:size(Z,2)
        num=size(Z,2)-zj+1;
        for ni=1:num
            zi=size(Z,1)-ni+1;
            tmp(ni)=Z(zi,zj+ni-1);
        end
        tmp=tmp(1:num);
        flat_z_pred(zj-1)=mean(tmp);
        sd_flat_z_pred=sd_flat_z_pred+std(tmp);
    end
    sd_flat_z_pred=sd_flat_z_pred/(size(Z,2)-1);

    hz=hankel(flat_z,[flat_z(m) flat_z_pred]);
    %% temporal component exaction ,  Y=U1S1V1, part 1
    [U1,S1,V1]=svd(hz');
    pc_dpca=V1';
    pc_num=3;

    variance=abs(diag(S1));
    variance_dpca=variance/(sum(variance)+0.0001);

    [pc,expvar]=eig(hz'*hz);
    [expvar,idx]=sort(real(diag(expvar)),'descend');
    pc_dpca=real(pc(:,idx));

    Zdpca=hz*pc_dpca;
    pc_dpca=Zdpca';

    variance=abs(diag(expvar));
    variance_dpca=variance/(sum(variance)+0.0001);

%% W_comb
 [weight_W_1,xi]=fpca_rebuild(H,n,L);
weight_W_1=weight_W_1/norm(weight_W_1,'fro');
Z1=X*weight_W_1;

    if ii==start
        weight_Z_1=norm(Z1,'fro');
    end
    
    ii
    clear flat_z_1 flat_z_pred_1
    sd_flat_z_1=0;
    sd_flat_z_pred_1=0;
    for zi=1:size(Z1,1)
        num=0;
        for zj=1:size(Z1,2)
            if zi-zj+1<1
                break;
            end
            num=num+1;
            tmp(num)=Z1(zi-zj+1,zj);
        end
        tmp=tmp(1:num);
        flat_z_1(zi)=mean(tmp);
        sd_flat_z_1=sd_flat_z_1+std(tmp);
    end
    sd_flat_z_1=sd_flat_z_1/size(Z1,1);
    
    for zj=2:size(Z1,2)
        num=size(Z1,2)-zj+1;
        for ni=1:num
            zi=size(Z1,1)-ni+1;
            tmp(ni)=Z1(zi,zj+ni-1);
        end
        tmp=tmp(1:num);
        flat_z_pred_1(zj-1)=mean(tmp);
        sd_flat_z_pred_1=sd_flat_z_pred_1+std(tmp);
    end
    sd_flat_z_pred_1=sd_flat_z_pred_1/(size(Z1,2)-1);

    hz_1=hankel(flat_z_1,[flat_z_1(m) flat_z_pred_1]);
    %% temporal component exaction ,  Y=U1S1V1, part 1
    [U1_1,S1_1,V1_1]=svd(hz_1');
    pc_dpca_1=V1_1';
    pc_num=3;
    
    variance_1=abs(diag(S1_1));
    variance_dpca_1=variance_1/(sum(variance_1)+0.0001);
 
    [pc_1,expvar_1]=eig(hz_1'*hz_1);
    [expvar_1,idx_1]=sort(real(diag(expvar_1)),'descend');
    pc_dpca_1=real(pc_1(:,idx_1));
    
    Zdpca_1=hz_1*pc_dpca_1;
    pc_dpca_1=Zdpca_1';
    
    variance_1=abs(diag(expvar_1));
    variance_dpca_1=variance_1/(sum(variance_1)+0.0001);


embedding_error = norm(Z-hz, 'fro')/norm(Z, 'fro')
embedding_error_1= norm(Z1-hz_1, 'fro')/norm(Z1, 'fro')
norm(Z(1:m-1,2:L)-Z(2:m,1:L-1),'fro')/norm(Z(1:m-1,2:L), 'fro')-norm(Z1(1:m-1,2:L)-Z1(2:m,1:L-1),'fro')/norm(Z1(1:m-1,2:L), 'fro')

embedding_error_about_X(1,ii-155)=embedding_error;
embedding_error_about_X(2,ii-155)=embedding_error_1;
 if embedding_error>embedding_error_1||embedding_error==embedding_error_1
    num_rebuild_nice1=num_rebuild_nice1+1;
end

end
%% mse
mean_embedding_error(1)=mean(embedding_error_about_X(1,:))
mean_embedding_error(2)=mean(embedding_error_about_X(2,:))


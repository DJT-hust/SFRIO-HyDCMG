function [omega, mu, rho]=Construcrt_Box_uncertainty_set(err,epsilon,delta)

data=err;

n_data=size(data,2);
N=size(data,1);
n2=ceil(log(delta)/log(1-epsilon));
n1=n_data-n2;
[data_set1,data_set2]=Seperate_dataset(data,n1,n2);

mu=mean(data_set1');
D=cov(data_set1');
Q=inv(D);
[P,LAMBDA]=eig(Q);
sqrt_LAMBDA=sqrt(LAMBDA);

projection_max=zeros(n2,1);
omega=sqrt_LAMBDA*P;
for k=1:n2
map_pos=omega*(data_set2(:,k)-mu');
all_map=[map_pos; -map_pos];
projection_max(k)=max(all_map);
end
projection_max_sort=sort(projection_max);
sum_p=0;
for k=0:n2
    sum_p=sum_p+nchoosek(n2,k)*(1-epsilon)^k*epsilon^(n2-k);
    if sum_p>1-delta
        k=k+1;
        break;
    end
end
rho=projection_max_sort(k);
omega=[omega;-omega];
end

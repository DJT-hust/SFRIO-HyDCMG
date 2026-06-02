 function [data_set1,data_set2]=Seperate_dataset(data,n1,n2)

n_data=n1+n2;
rowrank = randperm(n_data);
rowrank_data_set1=rowrank(:,1:n1);
rowrank_data_set2=rowrank(:,n1+1:n_data);
data_set1=data(:,rowrank_data_set1);
data_set2=data(:,rowrank_data_set2);

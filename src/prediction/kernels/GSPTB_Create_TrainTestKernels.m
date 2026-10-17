function [Ker_train,Ker_test] = GSPTB_Create_TrainTestKernels(Data,id_train,id_test)
    Ker_train = corr(Data(:,id_train));
    Ker_test = corr(Data(:,id_train),Data(:,id_test))';
end
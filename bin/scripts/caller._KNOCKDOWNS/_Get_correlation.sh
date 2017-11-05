#!/bin/bash



%The batch effect could be easily?? corrected by comparing the rank of log2 transformed transcripts
%The siGFP correlation increases from 0.9435 to 0.9813

file='Normalization/293wt-RNAseq_HEK293_Nanofectin_RNAseq_HEK293_siAUF1_RNAseq_HEK293_siGFP_RNAseq_HEK293_siTIA1_RNAseq_HEK293_sihnRNPC_RNAseq_mRNASeq_No_4SU_No_XL_rep_A_Clip13_mRNASeq_No_4SU_No_XL_rep_B_Clip13_si_GFP_mRNASEQ_si_HuR_mRNASEQ.raw';
file='without_ribosomal_proteins.raw';
p=10e-4
p=10e-6

a=importdata(file);
a2=a.data;
a2t=sum(a2);
al=log2(a2./repmat(a2t, size(a2,1),1)+p);
%al=log2(a2./repmat(a2t, size(a2,1),1)+p); %lower df degrees of freedom?
corr(al)
figure, imagesc(corr(al));

corr((a.data(:,4)), (a.data(:,10)))      
%    0.9435
corr(log2(a.data(:,4)+p), log2(a.data(:,10)+p), 'type', 'spearman')
%    0.9650
corr(log2(a.data(:,4)+p), log2(a.data(:,10)+p))
%    0.9527  %w p=10e-4
corr(al(:,1),al(:,2))
%    0.9675  %w p=10e-4
%    0.9813  %w p=10e-6

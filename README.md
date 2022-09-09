# Genome-wide inference reveals that feedback regulations constrain promoter-dependent transcriptional burst kinetics


Code and data for inference and analysis of feedback regulations type and burst kinetics from single-cell RNA-seq described in the article 'Genome-wide inference reveals that feedback regulations constrain promoter-dependent transcriptional burst kinetics' by Songhao Luo, Zihao Wang, Zhenquan Zhang, Tianshou Zhou, and Jiajun Zhang. biorxiv doi: https://doi.org/10.1101/2022.04.08.487618

## System Requirements

### **Scripts of inference**

The code of parameter inference for the hierarchical and telegraph models in the `code_Inference` folder are all MATLAB scripts with versions

```matlab
MATLAB Version: 9.10.0.1649659 (R2021a) Update 1
Operating System: Microsoft Windows 10 Professional Version 10.0 (Build 19044)
Java Version: Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
```

### **Scripts of data analysis and visualization**

The code of data analysis and visualization for investigating the burst kinetic in the `code_analysis&plot` folder are all R scripts with versions

```R
R version 4.1.3 (2022-03-10) -- "One Push-Up"
Copyright (C) 2022 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64 (64-bit)
```

## Directories

### **code_Inference**

#### Hierarchical_model

**master_demo.m**

*master_demo.m* is the master program for estimating the burst kinetic parameters and feedback type from synthetic data simulated by random sampling for a given parameter. You can set any parameters ($a,b,k,h$) to test the accuracy of the inference algorithm.

output : *master_demo.m* outputs the vector of estimated parameters, result plot of model selection. All the data generated during the inference process is stored by a mat file in `results\validation_example\` directory.

**master_MEFs.m**

*master_MEFs.m* is the master program for mouse embryonic fibroblasts (MEFs) scRNA-seq data. The program can be accelerated using parallel programs. The results for each individual gene are stored in `results/results_MEF` directory and csv.file of all results can be collected by *collectResult.m*. For example:

|      a      |      b      |     k      |           h            | ...  |
| :---------: | :---------: | :--------: | :--------------------: | ---- |
| BF of gene1 | BS of gene1 | k of gene1 | Feedback type of gene1 | ...  |
| BF of gene2 | BS of gene2 | k of gene2 | Feedback type of gene2 | ...  |
|     ...     |     ...     |    ...     |          ...           | ...  |

**validation_Feedback.m/validation_NonFeedback.m/validation_CellNumber.m/validation_Sensitive.m**

Some code about validating the efficiency of the inference algorithm

**Directory**

The `autoInference` directory contains scripts for inferring burst kinetic parameters and feedback type.

The `data` directory contains the scRNA-seq data from mouse embryonic fibroblasts described in Larsson  et al. 2019, Nature. MEF_QC_all is quality controlled and genename.txt is the gene name corresponding to each row.

The `results` directory contains results about demo example, validation, MEFs. 

The `misc` directory contains some other functions.

**Example**

```matlab
# open the test_demo.m
clear;clc;
addpath(genpath(pwd))
a = 2;
b = 10;
k = 10;
h = 3;
e = 0.05;
r = 0.5;
theta_true = [a,b,k,h];

% Generating simulated data
m_max = round(r*a*b + 10*sqrt(r*a*b+r^2*a*b^2));
[data,p_true] = generateSample(theta_true,m_max,50000);
data_mean = mean(data);
data_var = var(data);
data_noise = data_var/data_mean^2;

% Figure the data and p_true
figure
histogram(data,max(data),'normalization','pdf')
hold on
plot(0:length(p_true)-1,p_true)
hold off

%% Inference
% Inference for no feedback
results_non_total= [];
lb_non = [1e-1, 1];% Lower bound
ub_non = [30 , 20 ];% Upper bound
if data_mean < data_var
    results_non_total = inferenceKinetic(data,lb_non,ub_non,'non-feedback');
end

% Inference for positive feedback (H<0)
lb_positive = [1e-1, 1,    1, -10];% Lower bound
ub_positive = [30 , 20 ,  1e3,  -1];% Upper bound
results_positive_total = inferenceKinetic(data,lb_positive,ub_positive,'feedback');

% Inference for negative feedback (H>0)
lb_negative = [1e-1, 1,    1,   1];% Lower bound
ub_negative = [30 , 20 ,  1e3,  10];% Upper bound
results_negative_total = inferenceKinetic(data,lb_negative,ub_negative,'feedback');

% Model selection
results_total = [results_non_total;results_positive_total;results_negative_total];
theta_est = modelSelect(data,results_total,lb_non,ub_non);
save(sprintf('results\\validation_example\\validation_example_a=%d_b=%d_k=%d_h=%d.mat',a,b,k,h))
```

#### OnOff_model

 **master.m**

*master_MEFs.m* is the master program for mouse embryonic fibroblasts (MEFs) scRNA-seq data based on the telegraph model. The program can be accelerated using parallel programs. The results for each individual gene are stored in `results/results_MEF` directory and csv.file of all results can be collected by *collect_result.m*.

**Other scripts**

Other matlab files are the necessary scripts for inferring.

**Directory**

The `data` directory contains the scRNA-seq data from mouse embryonic fibroblasts described in Larsson  et al. 2019, Nature. MEF_QC_all is quality controlled and genename.txt is the gene name corresponding to each row.

The `results` directory contains results about demo example, MEFs. 

### **code_analysis&plot**

This folder contains a number of subfolders with corresponding folder names for the figures that appear in the manuscript and supplementary files. Each subfolder contains the analysis code and visualization code to generate the images as well as the pdf that has been generated. Before running the code, please open `genomeAnaly.Rproj` to ensure the accuracy of the path.


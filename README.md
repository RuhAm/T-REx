*T-REx*
Tensor decomposition-based Robust feature Extraction and classification to detect natural selection from genomic data.

*T-REx* is a machine learning algorithm that is able to differentiate between signatures of sweep from neutrality. We apply tensor decomposition to images of haplotypes across sampled individuals, and then use these extracted features as input to classical linear (elastic net logistic regression [EN]) and nonlinear (support vector machine with a radial basis kernel [SVM] and random forest [RF]) models. Finally, we find the predicted probabilities of sweep along with the predicted class for every test observation. In this software package, we have outlined the step by step procedure towards using *T-REx* in order to detect natural selection. In addition to that, we have included data preprocessing steps that can be used to parse simulated data that are in *.ms* format and empirical data in *.vcf* format.


## 1.0 Download and installation

To download and install *T-REx*, please follow the instructions below:

Download the repository using git.
```
git clone https://github.com/RuhAm/T-REx.git
```

*T-REx* uses R version 4.2.1 for the machine learning pipeline and Python3 for the data preprocessing. 



The following libraries need to be installed in *R* environment.
```
install.packages(c("abind", "MASS", "glmnet","rTensor","caret","ranger","tidyverse","ggpubr"))
```
```
install.packages("liquidSVM", repos="http://pnp.mathematik.uni-stuttgart.de/isa/steinwart/software/R")
```
The following libraries need to be installed in a *Python3* environment.
```
pip install pandas numpy math scipy skimage scipy argparse
```


## 2.0 Training and Testing from simulated data

*T-REx* uses data in *.ms* format and converts them into haplotypes which is further downsampled into 64*64 dimensional matrix using the following commands.

2.1 Once the repository is downloaded, please open a terminal at the *T-REx* directory. This can be acheived using the following commands given the repository is downloaded in the desktop folder. 


```bash
$ cd /Users/user/Desktop/T-REx
```

2.2 Please use the following commands to preprocess the data in *.ms* format for sweep and neutral observations to be used for training:

2.2.1 The "Data" folder contains a folder named "MS_files_train". The *.ms* files for sweep and neutral observations need to be placed in this folder. There are 100 example files located in the folder. The sweep files should are named as follows:

```bash
sweep_1.ms, sweep_2.ms ... sweep_100.ms
```

The neutral files should are named as follows:

```bash
neut_1.ms, neut_2.ms ... neut_100.ms
```


2.2.2 To preprocess the *.ms* files into *.csv* format please use the following commands.

   ```sh
   $ Python3 train_ms.py <number of files> <class>
   ```
The first argument <number of files> is the number of files that the user wants to preprocess. The second argument  <class> which takes on 0 or 1 as value. 1 is going to preprocess sweep observations and 0 is going to preprocess neutral observations.

2.2.3 To preprocess the *.csv* files using our unique alignment processing strategy, please use the following commands:

   ```sh
   $ Python3 parse_train.py <number of files> <class>
   ```
The first argument <number of files> is the number of files that the user wants to preprocess. The second argument is <class>, which takes on 0 or 1 as value. 1 is going to preprocess sweep observations and 0 is going to preprocess neutral observations.

The preprocessed files are going to be located at the following directory:

   ```sh
   "/Users/user/Desktop/T-REx/Data/CSV_files/"
   ```



2.3 Please use the following commands to preprocess the data in *.ms* format for sweep and neutral observations to be used for testing:


2.3.1 The "Data" folder contains a folder named "MS_files_test". The *.ms* files for sweep and neutral observations need to be placed in this folder using a prefix that are shown below. There are 100 example files located in the folder. The sweep files should are named as follows:

```bash
sweep_1.ms, sweep_2.ms ... sweep_100.ms
```

The neutral files should are named as follows:

```bash
neut_1.ms, neut_2.ms ... neut_100.ms
```


2.3.2 To preprocess the *.ms* files into *.csv* format please use the following commands.

   ```sh
   $ Python3 train_ms.py <number of files> <class>
   ```
The first argument <number of files> is the number of files that the user wants to preprocess. The second argument is the <class> which takes on 0 or 1 as value. 1 is going to preprocess sweep observations and 0 is going to preprocess neutral observations.

2.3.3 To preprocess the *.csv* files using our unique alignment processing strategy, please use the following commands:

   ```sh
   $ Python3 parse_test.py <number of files> <class>
   ```
The first argument <number of files> is the number of files that the user wants to preprocess. The second argument is the <class> which takes on 0 or 1 as value. 1 is going to preprocess sweep observations and 0 is going to preprocess neutral observations.



2.3.4 The preprocessed files are going to be located at the following directory:

   ```sh
   "/Users/user/Desktop/T-REx/Data/CSV_files/"
   ```


The train files would be saved in the following pattern:

```bash
sweep_train_align_1.csv, sweep_train_align_2.csv ...
neut_train_align_1.csv, neut_train_align_2.csv ...
```

The test files would be saved in the following pattern:
```bash
sweep_test_align_1.csv, sweep_test_align_2.csv ...
neut_test_align_1.csv, neut_test_align_2.csv ...
```


2.4 To perform tensor decomposition using a rank specified by the user, please use the following command:

```bash
$ Rscript TD.R <rank> <number of sweep train sample> <number of neutral train sample> <number of sweep test sample> <number of neutral test sample>
```

This command will perform tensor decomposition using a rank supplied by the user and train 3 classifiers (Elastic Net, Random Forest, Support Vector Machine) using example training data and output the probabilities of sweep using the example test data. The final output files will be saved in the following directory:

   ```sh
   "/Users/user/Desktop/T-REx/Data/Results/"
   ```



2.4.1 For Elastic net- the output probabilities will be saved in the "Probs_EN.csv" file and the class prediction will be saved in the "Class_EN.csv" file.

Similarly, output probabilities from Random Forest classifier will be saved in the "Probs_RF.csv" file and the class prediction will be saved in the "Class_RF.csv" file. 

Finally,  output probabilities from SVM will be saved in the "Probs_SVM.csv" file and the class prediction will be saved in the "Class_SVM.csv" file. 

## 3.0 Working with empirical data (VCF files)
3.1 To work with empirical data we need to keep the *.vcf* formatted file in the "VCF" folder inside the "Data" folder. One example file for chromosome 22 is given. The files should be named as follows: 

```bash
CEU(chromosome number).vcf
```
For example, 
```bash
CEU22.vcf
```




3.1.1 The *.vcf* file needs to be converted into *.ms* files using the following command.

   ```sh
   $ Python3 VCF_ms.py <name of the VCF file> 
   ```
For example, the following command will parse the "CEU22.vcf" file.

   ```sh
   $ Python3 VCF_ms.py CEU22
   ```
3.1.2 the *.ms* files will be saved in the following directory:

   ```sh
   "/Users/user/Desktop/T-REx/Data/VCF/MS_files/"
   ```
3.1.3 To preprocess the *.ms* files into *.csv* format please use the following commands:

   ```sh
   $ Python3 test_vcf.py
   ```

3.1.4 To preprocess the *.csv* files using our unique alignment processing strategy, please use the following commands:

   ```sh
   $ Python3 Parse_vcf.py <number of files>
   ```
The first argument <number of files> is the number of files that the user wants to preprocess.


The preprocessed files are going to be located at the following directory:

   ```sh
   "/Users/user/Desktop/T-REx/Data/VCF/CSV/"
   ```


3.2 To perform tensor decomposition using a rank specified by the user, please use the following command:

```bash
$ Rscript TD_vcf.R <rank> <number of sweep train sample> <number of neutral train sample> <number of test samples>
```




This command will perform tensor decomposition using a rank supplied by the user and train 3 classifiers (Elastic Net, Random Forest, Support Vector Machine) using example training data and output the probabilities of sweep using the empirical test data. 

The results will be saved in the following folder:

   ```sh
   "/Users/user/Desktop/T-REx/Data/Results/Empirical/"
   ```

3.2.1 For Elastic net- the output probabilities will be saved in the "Probs_EN.csv" file and the class prediction will be saved in the "Class_EN.csv" file.

Similarly, output probabilities from Random Forest classifier will be saved in the "Probs_RF.csv" file and the class prediction will be saved in the "Class_RF.csv" file. 

Finally,  output probabilities from SVM will be saved in the "Probs_SVM.csv" file and the class prediction will be saved in the "Class_SVM.csv" file. 




## 4.0 Contribution

This research work is supported by National Institutes of Health and National Science Foundation. 


Authors: Md Ruhul Amin, Mahmudul Hasan, Sandipan Paul Arnab, Michael DeGiorgio


## 5.0 Contact 

Md Ruhul Amin

aminm2021@fau.edu





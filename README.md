# Cone-OMP - Cone Orthogonal Matching Pursuit

We present a novel method used for anomaly detection. In the following scenario, we obtain good results in the detection of abnormal heartbeats. The method consists of the use of sparse representations which use atoms that are not vectors, but cones. We named this method Cone Orthogonal Matching Pursuit (Cone-OMP).

## Citing Us

This represents the main resources necessary for the reproduction of the results presented in “Sparse representations with cone atoms”, by Denis C. Ilie Ablachim, Andra Băltoiu, and Bogdan Dumitrescu. If you use our work please cite the following paper

```
@inproceedings{ilie2023sparse,
  title={Sparse Representations with Cone Atoms},
  author={Ilie-Ablachim, Denis C and Baltoiu, Andra and Dumitrescu, Bogdan},
  booktitle={ICASSP 2023-2023 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)},
  pages={1--5},
  year={2023},
  organization={IEEE}
}
```

## Project structure

```
"data" directory
	|_ dicts = directory containing the trained dictionaries, with K-SVD
	|_ preproc = directory containing preprocessed signals (on segments) from the MITDB dataset
"classif" - a general function which integrates all the classification methods
"classif_generate_roc" - function computing the proposed OMP method and the performance
"K_SVD" - DL method used for training dictionaries
"mexOMP" - function used for the standard OMP algorithm
"omp_cone_single" - function used for the Cone-OMP algorithm
"preproc_mitdb" - function used for preprocessing the signals available in the MTIDB dataset
"read_mitdb" - function used to download MITDB dataset
"test_methods" - main script for the test
```

## Reproducing the results

1. Download the Waveform Database Software Package (WFDB) for MATLAB and Octave. This is available [here](https://physionet.org/content/wfdb-matlab/0.10.0/). Download the zip, extract it and copy the resulting directory into the project. After this go in Matlab into the mcode directory available in the toolbox and run the following commands

        [old_path]=which('rdsamp'); if(~isempty(old_path)) rmpath(old_path(1:end-8)); end
        wfdb_url='https://physionet.org/physiotools/matlab/wfdb-app-matlab/wfdb-app-toolbox-0-10-0.zip';
        [filestr,status] = urlwrite(wfdb_url,'wfdb-app-toolbox-0-10-0.zip');
        unzip('wfdb-app-toolbox-0-10-0.zip');
        cd mcode
        addpath(pwd)
        savepath

2. Run the **read_mitbih** script to download the ECG signal from MIT-BIH Arrhythmia Database. By default, only record #109 will be downloaded but you can uncomment the necessary line to download the whole dataset.

3. Run **preproc_mitdb(109)** to preprocess the signals. By processing, we divide each signal into segments of 5 minutes that will be processed independently. The dimension of each segment is reduced from 256 to 32 by projecting them on the 32 leading PCA basis vectors of the segment.

4. Run a test by using the **test_methods** script. This will train a dictionary by the K-SVD algorithm, which implies the OMP method for sparse representation. After the training is complete, the classification stage starts. If you want to add new sparse representation methods you should edit the "classif" function. By default, we use the OMP and Cone-OMP methods.
We compute the ROC AUC and the number of FP (False Positive) results at true positive rate 100%, 97%, and 94%.
For record #109, this corresponds to TP=40, 39, and 38, respectively.


## Funding

This work is supported by a grant of the Ministry of Research, Innovation and Digitization, CNCS – UEFISCDI, project number PN-III-P4-PCE-2021-0154, within PNCDI III.

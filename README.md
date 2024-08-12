# BMAS Beamforming for plane wave imaging
Beam Multiply and Sum: A Tradeoff between Delay and Sum and Delay Multiply and Sum Beamforming for Ultrafast Ultrasound Imaging

Repository to share the MATLAB implementation of beam multiply and sum beamforming for plane wave ultrasound imaging.

###### DATE 31-07-2022 : VERSION 1.0
###### Last Updated on: Aug 2024

#### AUTHORS: 
A. N. Madhavanunni and Mahesh Raveendranatha Panicker <br />
Center for Computational Imaging, Indian Institute of Technology Palakkad, Kerala, India.

#### Instructions for the execution of codes 

1. Download the source_code and data folder from the repository.
2. Please note that the source code has the following dependencies:
   (a) PICMUS evaluation scripts: apodization.m, us_dataset.m and band_pass.m
      Download the codes from https://www.creatis.insa-lyon.fr/Challenge/IEEE_IUS_2016/sites/www.creatis.insa-lyon.fr.Challenge.IEEE_IUS_2016/files/archive_to_download.zip and place the above-mentioned files in source_code/lib/ before running the codes.

   (b) PICMUS dataset: carotid_cross_expe_dataset_rf.hdf5
      Download the dataset from https://www.creatis.insa-lyon.fr/Challenge/IEEE_IUS_2016/sites/www.creatis.insa-lyon.fr.Challenge.IEEE_IUS_2016/files/in_vivo.zip and place 'carotid_cross_expe_dataset_rf.hdf5' in data/ before running the codes.



#### Instructions for the execution of MATLAB App 

Please follow the below instructions to run the code.
1. Download the MATLAB app from the following link: https://drive.google.com/drive/folders/16GTG-9JxCSyHSQWgUPHUm_QGFJDvMIvV?usp=sharing
2. Install the app using BMAS_AppInstaller_web.exe in the package.
3. Run the BMAS_for_PWI.exe 

This application allows you to:
      (i) Vary the hyper-parameters (L and s) and reconstruct the B-mode images for various datasets using the BMAS beamformer.
      (ii) Tune the hyperparameters and see a transition from DAS beamforming to DMAS beamforming for the selected datasets.


##### Code Available under Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0) (see https://creativecommons.org/licenses/by-nc-nd/4.0/)

##### Please see the articles mentioned in Academic references and Acknowledgements that are to be cited for any usage of the code and/or data.

#### ACADEMIC REFERENCES TO BE CITED:
Details of the BMAS beamforming is available in the following articles by A. N. Madhavanunni and Mahesh Raveendranatha Panicker <br />
1. A. N. Madhavanunni and Mahesh Raveendranatha Panicker. [Beam multiply and sum: A tradeoff between delay and sum and filtered delay multiply and sum beamforming for ultrafast ultrasound imaging](https://doi.org/10.1016/j.bspc.2023.104807). _Biomedical Signal Processing and Control_ 85 (2023): 104807.<br />

2. A. N. Madhavanunni and Mahesh Raveendranatha Panicker, [Lesion Detectability and Contrast Enhancement with Beam Multiply and Sum Beamforming for Non-Steered Plane Wave Ultrasound Imaging](https://doi.org/10.1109/ISBI52829.2022.9761662). _IEEE 19th International Symposium on Biomedical Imaging (ISBI)_, Kolkata, India, 2022.

3. A. N. Madhavanunni and Mahesh Raveendranatha Panicker, [Performance Evaluation of Beam Multiply and Sum Beamforming with Coherent Plane Wave Compounding: In-vitro Results](https://doi.org/10.1109/SAUS61785.2024.10563536). _IEEE South Asian Ultrasonics Symposium 2024 (SAUS 2024)_, Gandhinagar, India. March. 2024. 

#### Acknowledgements:
##### Ultrasound dataset: Plane-wave imaging challenge for medical ultrasound (PICMUS) of 2016 IEEE International Ultrasonics Symposium 
1. H. Liebgott, A. Rodriguez-Molares, F. Cervenansky, J. A. Jensen, and O. Bernard, “Plane-Wave Imaging Challenge in Medical Ultrasound,” _IEEE International Ultrasonics Symposium, IUS_, vol. 2016-November, 2016.

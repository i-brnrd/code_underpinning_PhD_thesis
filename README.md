# code_underpinning_PhD_thesis
This repository contains the code underpinning Isla Barnard's PhD Thesis:

'Computational simulations of ultraviolet radiation penetration into human skin', 2021, University of St Andrews. 

![alt text](https://github.com/i-brnrd/code_underpinning_PhD_thesis/blob/main/images_mainrepo/smaller.png?raw=true)

## Introduction

This repository was created to allow creation of a single DOI for the code underpinning the thesis. As such, this single repo contains code from several different individual repositories- however, in order to prevent conflicts, the git repos themselves have been deleted. Each repo, as it stood at publication, resides in a directory in this repository. Code found within this repository has not been maintained since December 2020 or earlier. For more up to date, maintained code, please see publicly available repos by i-brnrd. 

Data produced by this code has been deposited in the PURE repository, found here: 
https://doi.org/10.17630/b379c3f3-cc1c-497b-bd28-6c7dc0cc209d 

This repo was created prior to publication of the thesis and as such, this README does not contain a direct link to full text of the thesis.
To source the full text, please search the St Andrews Research Repository: https://research-repository.st-andrews.ac.uk/
for the thesis titled, 'Computational simulations of ultraviolet radiation penetration into human skin' (also search for 'Monte Carlo' for other publications from our group!). 
If this fails please raise an issue here, and I should be notified  (I'm also available on researchgate: https://www.researchgate.net/profile/Isla-Barnard and on linkedIn). 
	 
## Code
All computational code is in Fortran 90; compiled using gfortran, run on an 8-core Intel® Xeon(R) CPU E3-1270 v5 @ 3.60GHz. 

All code is based on Dr. Kenny Wood's F77 MCRT grid code, found here: http://www-star.st-and.ac.uk/~kw25/teaching/mcrt/mcrt.html

### Chapters 2, 3 & 4 : cpd_thesis_dna
The code described in **Chapters 2 & 3** of the thesis is available, in some form, in any of the directories in this repository- however the code that is least altered from that as described in **Chapters 2 & 3**, and the code used in **Chapter 4**, is found in **cpd_thesis_dna**. Any users should start here as this is the UV-MCRT skin model in it's simplest form.


### Chapters 5 & 6: melanin_sunscreen
The alterations to the code as described in **Chapters 5 & 6** are found in **melanin_sunscreen**.


### Chapter 7: psoralen
The code as described in **Chapter 7** is found in **psoralen**.

### Chapter 8: UVC
The code as described in **Chapter 8** is found in **UVC**.


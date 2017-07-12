
### Picking a list of transitions for SRM

With my list of interesting proteins, I hit Skyline with the goal of compiling a list of ~16 proteins, each with 9 transitions that describe 3 peptides (3 transitions per peptide), which are consistently detectable in all samples, and the protein peaks look good. 

Tips from Yaamini: 
  * Delete the precursor ions from my list prior to assessing quality, to remove visual noise
  * Run through my entire list of ~200 proteins, since very few will meet the QC standards
  * Save-as Skyline project, then delete proteins until I have a skyline doc. with only my proteins/peptides/transitions of interest for SRM, then export. 

Okay, I have to start somewhere. 

First protein checked: comp134820_c0_seq2 | Pre-mRNA-splicing ATP-dependent RNA helicase PRP28 (EC 3.6.4.13)
Only has 1 peptide, with only 2 good transitions. Not consistent in all samples.
![image](https://user-images.githubusercontent.com/17264765/28047722-f2846fca-65a0-11e7-9335-46cc1eb1ff21.png)

Gigasin, unfortunately, looks like crud: 
![image](https://user-images.githubusercontent.com/17264765/28049188-7b9da728-65aa-11e7-9166-fa10b5800704.png)

Heat Shock Protein HSP 90-alpha | comp132209_c0_seq1: excellent
![image](https://user-images.githubusercontent.com/17264765/28049433-4a681bf0-65ac-11e7-9157-9320bbd6f4c5.png)
![image](https://user-images.githubusercontent.com/17264765/28049463-6742153c-65ac-11e7-8e26-f79c80f2d6c5.png)

Beta-1,3-glucan-binding protein (GBP) | comp135114_c0_seq1: poor
![image](https://user-images.githubusercontent.com/17264765/28049529-dd87660c-65ac-11e7-9016-78ba806c76df.png)

Arachidonate 5-lipoxygenase (5-LO) (5-lipoxygenase) (EC 1.13.11.34) | comp135856_c0_seq2: excellent
![image](https://user-images.githubusercontent.com/17264765/28049626-6c432458-65ad-11e7-99b4-c5a782e2a676.png)

I reviewed ~160 proteins from three lists: diff. expressed, stress-related, and Yaamini's targets. I counted the # peptides detected and assined a quality ranking: poor, fair, okay, good, great, excellent. Total # proteins with >2 peptides and ranked great or excellent: 40

![image](https://user-images.githubusercontent.com/17264765/28093795-cfe9d6c8-664d-11e7-8319-b3cae9cfd1e9.png)

From this list I selected the following: 
  * TBD

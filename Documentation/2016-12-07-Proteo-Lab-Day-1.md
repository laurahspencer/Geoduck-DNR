---
layout: post
title: Protein Extraction: day 1
---

## DNR Geoduck OA overview

Gearing up for the DNR Geoduck Ocean Acidification project's first, exploratory protein extractions.

Here is a short summary of the DNR Geoduck OA project:

The Pacific geoduck (Panopea generosa) is the largest clam native to the Pacific Northwest and is a burgeoning aquaculture species due to growing export demands from Asia. In Washington State, geoduck support important commercial fisheries via farmed and wild populations in Puget Sound. As a sedentary, calcifying bivalve occupying mostly subtidal sediment, geoduck are likely to be impacted by climate stressors, which have already been documented as trending towards warmer, more acidic marine conditions. In summer 2016 we performed two trials in which sibling juvenile geoduck were out-planted in five sites throughout Washington State, each site containing cohorts placed inside and outside eelgrass beds. Geoduck were enclosed to minimize predation, water chemistry was continuously monitored, and after four weeks of exposure geoduck gill tissue was taken for proteomic analysis. 

The gill proteome will be sequenced using tandem mass spectrometry to reveal expression pattern associated with local conditions. With these results we seek to examine how protein profiles might provide valuable information on local conditions including how environmental change can influence bivalve physiology. 

For details more details, data, protocols, scripts see the [Github repo](https://github.com/RobertsLab/project-geoduck-oa)

---

## Lab day 1 tasks

I am following the [Protein Preparation for MSMS protocol](https://github.com/sr320/LabDocs/blob/master/protocols/ProteinprepforMSMS.md), provided by [Emma Timmins-Schiffman](https://github.com/emmats), and updated by [Rhonda Elliott](https://github.com/Ellior2). 

[Rhonda's lab notebook](https://ellior2.github.io/), particularly posts from 11/16/2016 - 12/6/2016 (e.g. ["First Lab Day"](https://ellior2.github.io/FirstLabDay/) from her recent extraction, is also very helpful.  
 
### Today's Tasks:
* Select samples for this run: 
* Prepare snap top centrifuge tubes - 3 per sample
* Get dry ice from biochem building  
* Cut samples in half to save one half for DNA analysis 
* Make ammonium bicarbonate

#### Select Samples:
* This first Mass Spec run is designated as "exploratory". 
* I am focusing on samples from Round #1, which was the first cohort of  geoduck outplanted for ~1 month. 
* There is only enough space in the Mass Spec to run 10 samples & 1 blank (the run is scheduled for late January 2017). Therefore, to get a broad swath, I opted for 2 samples per site (5 sites total), one from each the bare and eelgrass patches. This should provide a snapshot of every treatment in Round #1.
* To randomly select the sample numbers from each treatment I wrote an [R script](https://github.com/RobertsLab/project-geoduck-oa/blob/master/Trial1_docs/RandomSample.R). Samples selected are the following:  
  * G10 
  * G18
  * G48
  * G58
  * G68
  * G77
  * G92
  * G97
  * G119
  * G131 
  
#### Prepare snap top centrifuge tubes  
* In preparation for the sonication step tomorrow, and dividing my samples in 2 later today, I labeled 3 centrifuge tubes:
  * [Sample ID] A = 1/2 of raw sample pulled for proteomic Analysis
  * [Sample ID] B = BCA Assay, aka protein quantification
  * [Sample ID] C = Supernatant, which will become the sample for the Trypsin digestion step
  
#### Get dry ice from biochem building  
* Took ~20 minutes. Need: budget #, cooler.

#### Cut samples in half  
* Kept samples cold while processing by burying in dry ice
* Cleaned instruments via 10% bleach solution, rinsed with nano-pure water
* "Primed" weigh boats & new vial by placing them on the dry ice for ~1 minute, to cool
* Used tweezers & mini-scoops (whicever worked best, depending on the sample) to remove tissue from tube
* Using new razor blade for each sample, cut tissue in half. Returned half to the original sample tube, and the other to newly labeled "A" tubes.
* Returned original sample tube, now with half the sample, to the original, organized sample box in -80.

#### Prepare Ammonium Bicarbonate
* **This must be made for use within 24hrs**
* Made 20ml of the 50mM NH4HCO3 in 6M Urea:
  * Added 10ml of nanopure using a glass graduated cylinder into a 45ml falcon tube.
  * Weighed  79.06mg of NH4HCO3 in weigh boat, added to the falcon tube and vortexed to dissolve. 
  * Weighed  7.21g of Urea and added this to the falcon tube. Vortexed to dissolve. 
  * Poured the solution back into the graduated cylinder, then topped off volume with nanopure to total of 20ml. Poured solution back into falcon tube. Vortexted to dissolve.


 


  
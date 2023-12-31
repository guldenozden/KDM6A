 # Details of our experiments listed below;

**Cell Culture**
The BdEC (primary bladder epithelial cells) (ATCC® PCS-420-010™), SV-HUC-1 (immortalized uroepithelial cells) (ATCC® CRL-9520™), T24 (urinary bladder transitional 
cell carcinoma) (ATCC® HTB-4™) and HEK293T cell lines were used  for the experiments. The BdEC cell line was grown in the appropriate medium (ATCC® PCS-440-030™) 
and supplemented with the Prostate Epithelial Cell Growth Kit (ATCC PCS-440-040), as recommended by the ATCC, and Primocin (Invivogen). SV-HUC-1 was grown with 
RPMI-1640, and T24 and HEK293T were grown with DMEM, supplemented with 10% FBS and 1% penicillin. The cells were maintained in a humidified incubator at 37 ◦C and 5% CO2.

**Chromatin Immunoprecipitation**
The chromatin immunoprecipitation assays were performed by following the protocol of Weber et al., 2007. Briefly, the cells were fixed with 1% formaldehyde for 15 min
when they reached 70–80% confluence. A final concentration of 0.125 M glycine was used to prevent crosslinking. After the pellet was lysed with lysis buffer, the 
chromatin was fragmented using a S220 Covaris Ultrasound Sonicator. With the intent of determining the optimum sonication for 200–500 bp fragmented DNA, a 50 µL input sample
was taken from the lysate. The rest of the chromatin was used for the immunoprecipitation (IP). For IP, the chromatin was first pre-cleared using Dynabeads (Invitrogen, 
Waltham, MA, USA, 11203D) blocked with BSA and tRNA, and then, bound with 5 µL KDM6A antibody (Cell Signaling Technology, Danvers, MA, USA, 33510) overnight at 4 ◦C. 
The next day, the chromatin–antibody complex was bound with pre-blocked Dynabeads for 3 h at 4 ◦C, and after the washing steps, the chromatin was eluted. For both the 
input and IP DNA, following the RNase and proteinase K treatment, the samples were incubated overnight at 65 ◦C for reversal of the crosslink. The DNA was eluted using a
DNA Clean & Concentrator kit (Zymo Research, Orange, CA, USA, D4034).

**ChIP-Seq**
For the T24, BdEC, and SV-HUC-1 cell lines, ChIP and input DNA were sent to the EMBL GeneCore facility for library preparation and sequencing. A NEBNext DNA Ultra II kit 
was used for library preparation. The libraries were sequenced on a Nextseq 500 platform using 75bp SE high-output mode.

![image](https://github.com/guldenozden/KDM6A/assets/147516133/06ba0f0b-dca9-46a2-9610-f2e52ca64a5b)

**RT-qPCR**
Total RNAs were isolated using an MN Nucleospin RNA kit from each cell line. A total of 0.25 µg RNA was used for cDNA conversion using a Maxima First Strand cDNA Synthesis 
Kit. PCR was performed using the Applied Biosystems 7500 Fast Real-Time machine with a Roche FastStart Essential DNA Green Master (SYBR, Roche, Penzberg, Germany) kit for 
the selected genes. The experiments were performed as three technical replicates.

**ChIP-qPCR**
After performing the KDM6A ChIP experiments for the T24, SV-HUC-1, and BdEC cell lines, qPCR was performed using the Roche FastStart Essential DNA Green Master (SYBR) kit on 
an Applied Biosystems 7500 Fast Real-Time machine, and the results are presented as % input IP. The experiments were performed as three technical replicates. The ChIP-qPCR 
results were calculated using the ∆∆ct method. ChIP Ct values were normalized to the input.

**Primer Design**
The primers for RT-qPCR and ChIP-qPCR were designed using the primer designing tool Primer-Blast (NCBI). 

**Western Blotting**
The cells were washed with ice-cold PBS and scraped with 1× RIPA buffer supplemented with a protease inhibitor cocktail. After vortexing and incubating for 10 min on ice three 
times, the lysates were centrifuged at max speed at 4 ◦C for 20 min and the supernatants were collected. A Pierce BCA Protein Assay Kit (Thermo Fisher Scientific, Rockford, IL,
USA, 23225) was used to determine the protein concentrations. A total of 30 µg proteins were loaded onto 8% acrylamide:bis-acrylamide gel. The gel was run at 90 V for 30 min and 
at 120 V for one and a half hours. Gel transfer was performed at 350 mA for one and a half hours in a box full of ice. Non-fat dried milk powder was used for blocking the membrane 
for one hour at room temperature. KDM6A (CST 33510), HES1 (CST 11988), TLE1/2/3/4 (CST 4681S), and b-actin (CST 3700) primary antibodies were used with a 1/1000 dilution, and CST
5151S and LICOR 926-68072 secondary antibodies were used with a 1/30,000 dilution. Nitrocellulose membrane was incubated overnight at 4 ◦C with the KDM6A primary antibody, and for
one and a half hours at room temperature with the b-actin antibody. Secondary antibody incubations were performed for one hour at room temperature. The membrane was imaged using 
the Li-COR ODYSSEY Clx machine (LiCOR, Lincoln, NE, USA) in auto mode.

![image](https://github.com/guldenozden/KDM6A/assets/147516133/bf907e22-2573-4575-b9a5-a6f017b70659)

**Identification of KDM6A Mutation**
DNA was extracted from the cell lines BdEC, T24, and SV-HUC-1 using a ZymoResearch (Irvine, CA, USA) Quick DNA Miniprep Plus kit (D4068). Mutation analysis was performed using the 
Archer VariantPlex Myeloid Panel (Diagnostica Longwood, Zaragoza, Spain), an NGS panel containing the KDM6A gene. Mutation calling was performed  using ArcherDx Analysis (Version 6.2.7)
with default settings using the Gen-ERA NGS service.

![image](https://github.com/guldenozden/KDM6A/assets/147516133/3e228708-1f5d-4bbf-a20a-d981e258ebae)

**Plasmid Constructs and Transfection**
To generate the FLAG-KDM6A construct, a pCMV_HA_UTX (Addgene, Plasmid #24168, Watertown, MA, USA) vector was digested for 2 h at 37 ◦C using Kpn1 and Not1 restriction 
enzymes, and the KDM6A gene was removed from this vector. In order to add a FLAG-tag to KDM6A, 2-step PCR was set up using the Q5 polymerase (M0494S, Neb, Ipswich, MA, 
USA) with designed primers (see primer list, Table S3). The PCR product was cut from the agarose gel and isolated using NucleoSpin® Gel and a PCR Clean-up kit (MN, 740609, 
Düren, Germany). The digested FLAG-tag KDM6A gene and the digested pcDNA3.1(+)/myc-HisA (Invitrogen) vector were ligated for 2 h at room temperature using T4 Ligase (Neb, 
M0202S) in accordance with the Neb Ligation Calculator (https://nebiocalculator.neb.com/#!/ligation, accessed on 1 April 2022). The ligation product (pcDNA3.1_FLAG_KDM6A) 
was transformed into an E. coli 10beta bacterial strain. Growing colonies were isolated using the plasmid isolation kit (MN, 740588) and validation was performed via Sanger 
sequencing. A total of 80% confluent 10cm plate HEK293T cells were transfected with 7.5 µg plasmid DNA using Lipofectamine 3000 (Invitrogen, L3000015) reagent. The medium 
was changed 24 h after transfection. To generate the KDM6A truncated mutation like in the T24 cell line, a pcdna3.1_FLAG-KDM6A plasmid was digested for 2 h at 37 ◦C using
Kpn1 and Not1 restriction enzymes. After being run on 1% agarose gel at 100 V for 45 min, the FLAG-KDM6A band was isolated using NucleoSpin Gel and the PCR Clean-up kit 
(MN, Düren, Germany) from the gel. two-step PCR was set up using the Q5 polymerase with designed primers for E895* mutation. After this, the same steps as for the generation
of the pcDNA3.1_FLAG_KDM6A plasmid were followed.

**Co-IP**
Proteins from untransfected and pcDNA3.1_FLAG_KDM6A or pcDNA3.1_FLAG_E895*mutant_KDM6A transfected HEK293T cells (48 h post transfection) were isolated using 300 µL lysis buffer
(50 mM Tris-HCl pH:7.4, 150 mM NaCl, 1 mM EDTA, 1% TritonX100, Protease Inhibitor Cocktail (Roche, Indianapolis, IN, USA)). For input controls 20 µL untransfected and transfected 
HEK293T cell lysates were used. Remaining proteins were bound overnight at 4 ◦C on the rotator with Flag affinity gel (Sigma, A2220, St. Louis, MO, USA) in accordance with the 
manufacturer’s instructions. After the incubation, IP proteins were centrifuged, and pellets were washed with TBS three times. Proteins for both input and IP fractions were boiled 
with SDS Sample Buffer (62.5 nM Tris-HCl pH:6.8, 2%SDS, 10%(v/v) Glycerol and 0.002% Bromophenol Blue) at 95 ◦C for 3 min and analyzed by Western blotting. KDM6A (CST 33510), 
HES1 (CST 11988) and TLE1/2/3/4 (CST 4681S) antibodies were used for Western blotting.

![image](https://github.com/guldenozden/KDM6A/assets/147516133/484d2a74-abc4-4db7-b9bc-2a12deed6861)

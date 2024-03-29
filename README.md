
# **SNPMStat**

# **SNPMStat v4.0 : Statistical Analysis of SNP-Disease Association with Missing Genotype Data**

**SNPMStat** is a command-line program for the statistical analysis of
SNP-disease association in case-control/cohort/cross-sectional studies
with potentially missing genotype data. **SNPMStat** allows the user to
estimate or test SNP effects and SNP-environment interactions by
maximizing the (observed-data) likelihood that properly accounts for
phase uncertainty, study design and gene-environment dependence. For
SNPs without missing data, the program performs the standard association
analysis. For typed SNPs with missing data or untyped SNPs, the program
performs the maximum-likelihood analysis described in [Lin, Hu and Huang
(2008)](http://dlin.web.unc.edu/wp-content/uploads/sites/1568/2013/01/LinHuHuang08.pdf)
and Hu, Lin and Zeng (2010). We are working intensely to improve the
capabilities of **SNPMStat**, so please check back frequently for
updates.

## **SYNOPSIS**

**SNPMStat** \[**-sfile** specfile\] \[**-pfile** phenofile\]
\[**-gfile** genofile\] \[**-efile** exterfile\] \[**-h**\]
\[**-hfile** haplofile\] \[**-ofile** outfile\] \[**-no_fix**\]
\[**-no_remove**\] \[**-speed**\] \[**-ne**\] \[**-tag**\]
\[**-window**\] \[**-max**\] \[**-min\]**

## **OPTIONS**

 

<table style="width:99%;">
<colgroup>
<col style="width: 8%" />
<col style="width: 6%" />
<col style="width: 17%" />
<col style="width: 66%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Option</th>
<th style="text-align: left;">Parameter</th>
<th style="text-align: left;">Default</th>
<th style="text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;"><strong>-sfile</strong></td>
<td style="text-align: left;">{specfile}</td>
<td style="text-align: left;">specification.txt</td>
<td style="text-align: left;">Specify the specification file</td>
</tr>
<tr class="even">
<td style="text-align: left;"><strong>-pfile</strong></td>
<td style="text-align: left;">{phenofile}</td>
<td style="text-align: left;">phenotype.dat</td>
<td style="text-align: left;">Specify the phenotype file</td>
</tr>
<tr class="odd">
<td style="text-align: left;"><strong>-gfile</strong></td>
<td style="text-align: left;">{genofile}</td>
<td style="text-align: left;">genotype.dat</td>
<td style="text-align: left;">Specify the genotype file</td>
</tr>
<tr class="even">
<td style="text-align: left;"><strong>-efile</strong></td>
<td style="text-align: left;">{exterfile}</td>
<td style="text-align: left;">external.dat</td>
<td style="text-align: left;">Specify the external file</td>
</tr>
<tr class="odd">
<td style="text-align: left;"><strong>-h</strong></td>
<td style="text-align: left;"></td>
<td style="text-align: left;">No haplotype file</td>
<td style="text-align: left;">Use the haplotype file</td>
</tr>
<tr class="even">
<td style="text-align: left;"><strong>-hfile</strong></td>
<td style="text-align: left;">{haplofile}</td>
<td style="text-align: left;">haplotype.dat</td>
<td style="text-align: left;">Specify the haplotype file</td>
</tr>
<tr class="odd">
<td style="text-align: left;"><strong>-ofile</strong></td>
<td style="text-align: left;">{outfile}</td>
<td style="text-align: left;">results</td>
<td style="text-align: left;">Specify the output file. Unlike other file
names above, this one should be rid of any suffix.</td>
</tr>
<tr class="even">
<td style="text-align: left;"><strong>-no_fix</strong></td>
<td style="text-align: left;"></td>
<td style="text-align: left;">Perform internal check</td>
<td style="text-align: left;">Turn off internal check</td>
</tr>
<tr class="odd">
<td style="text-align: left;"><strong>-no_remove</strong></td>
<td style="text-align: left;"></td>
<td style="text-align: left;">Remove SNPs that cannot be aligned</td>
<td style="text-align: left;">Turn off default removal of SNPs that
cannot be aligned</td>
</tr>
<tr class="even">
<td style="text-align: left;"><strong>-speed</strong></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"><ul>
<li><p>1) Skip analysis of untyped SNPs that are not in the haplotype
file.</p></li>
<li><p>2) Exclude SNPs not in the haplotype file as predictors.</p></li>
<li><p>This will significantly speed up the program, but may lose
important untyped SNPs. This flag is meaningful only when
<strong>-h</strong> is used.</p></li>
</ul></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><strong>-ne</strong></td>
<td style="text-align: left;"></td>
<td style="text-align: left;">Use external panel</td>
<td style="text-align: left;">No external file</td>
</tr>
<tr class="even">
<td style="text-align: left;"><strong>-tag</strong></td>
<td style="text-align: left;">{#tags}</td>
<td style="text-align: left;">4</td>
<td style="text-align: left;">Specify the number of tag SNPs used to
impute the SNP of interest.</td>
</tr>
<tr class="odd">
<td style="text-align: left;"><strong>-window</strong></td>
<td style="text-align: left;">{win size}</td>
<td style="text-align: left;">50,000 (bp)</td>
<td style="text-align: left;">Specify the maximum distance to the
untyped SNP within which typed SNPs are identified as candidate
tags.</td>
</tr>
<tr class="even">
<td style="text-align: left;"><strong>-max</strong></td>
<td style="text-align: left;">{#SNPs}</td>
<td style="text-align: left;">20</td>
<td style="text-align: left;">Specify the maximum number of typed SNPs
identified as candidate tags</td>
</tr>
<tr class="odd">
<td style="text-align: left;"><strong>-min</strong></td>
<td style="text-align: left;">{#SNPs}</td>
<td style="text-align: left;">8</td>
<td style="text-align: left;">Specify the minimum number of typed SNPs
identified as candidate tags</td>
</tr>
</tbody>
</table>

 

By default, SNPMStat analyzes both typed SNPs in a study with
potentially missing data and untyped SNPs that are on an external panel.
For each SNP of interest, SNPMStat first identifies candidate tags
within a distance (**-window**) and then finds the predefined number of
SNPs that yields the largest *MD* measure (Nicolae, 2006, Genetic
Epidemiology, 30, 703-717). If the number of candidate tags within that
distance is less than the minimum (-min), the distance is enlarged until
the minimum number is met. If the number of candidate tags within the
distance exceeds the maximum (**-max**), only the closest maximum number
of SNPs are considered as candidate tags. We perform an internal check
to see whether the strand alignment between the study and external panel
can be determined from (a) the allele labels (at non A/T and G/C SNPs),
and (b) allele frequencies (at A/T and G/C SNPs). SNPs that cannot be
aligned are removed from the data. The internal check and removal can be
turned off by using **-no_fix** and **-no_remove**. Phased haplotypes
can be supplied by **-h** to facilitate the selection of tag SNPs. The
**-speed** flag can further speed up the selection process at the cost
of skipping SNPs that are not phased. If only typed SNPs are interested,
the use of external panel can be suppressed by **-ne**.

## **INPUT FILES**

### **specification file**

 

|                                 |
|---------------------------------|
| Example of a specification file |
| DESIGN = cohort                 |
|                                 |
| CATEGORICAL = smk_status        |
| DEPENDENT = smk_status CPD      |
|                                 |
| PANEL = 30 0 0                  |
|                                 |
| MODE = additive                 |
| EFFECT = G CPD G\*CPD           |
| OUTPUT = G G\*CPD               |

 

 

The specification file describes the feathers of the study and
variables, and specifies the disease risk model required for the
analysis. The syntax follows

> KEYWORD = value1 \[value2 …\]

with spaces around “=”. KEYWORD with an empty value, i.e., “KEYWORD =”,
is not allowed.

 

`DESIGN =`case-control/cohort/cross-sectional

> Specify the study design. **Required at the first line of the
> specification file.**

`CATEGORICAL =`{covariate names in the phenotype file}

> Specify covariates that are categorical (more than two levels). A
> categorical covariate is transformed into (level-1) indicators with
> the lowest level as the reference. For example, if *smoke* has values
> 1, 2, 3, it will be transformed into two indicators I(*smoke*=2) and
> I(*smoke*=3) with names “`smoke(2)`” and “`smoke(3)`“. Unspecified
> covariates are assumed to be continuous by default. **Optional.**

`DEPENDENT =`{covariate names in the phenotype file}

> Specify covariates that are potentially correlated with haplotypes.
> Unspecified covariates are assumed to be independent of haplotypes by
> default. **Optional.**

`PANEL =`  
{#trios \#duos \#singletons}

> Specify the number of trios, duos and singletons, respectively.  
> **  
> Required.  
> **

`MODE =`  
additive/recessive/dominant/codominant

> Specify the mode of inheritance. Default is additive mode.
> **Optional.**

`EFFECT =`{main effects and interactions}

> Specify the main effects and interactions considered in the disease
> risk model. In particular, the SNP effect is designated by ‘G’.
> Interactions between SNP and covariates are indicated by ‘\*’ with no
> space on either side. **Required.**

`OUTPUT =`{main effects and interactions}

> Specify the main effects and interactions whose estimation and testing
> results are to be outputted. These effects should be a subset of those
> in `EFFECT`. Each effect is outputted to a separate file, with the
> file name specified in **-ofile** appended by “`_effectname.out`“.
> Note that the “`*`” sign in interactions is replaced by “`$`” for
> legitimacy purpose. Specifying a categorical covariate induces
> multiple files corresponding to its derivative indicators. Specifying
> a codominant SNP effect induces two files corresponding to two
> genotypes. **Required.**

 

### **phenotype file**

 

|     |       |       |              |                             |
|-----|-------|-------|--------------|-----------------------------|
|     |       |       |              | Example of a phenotype file |
| Y   |   del |   age |   smk_status |   CPD                       |
| 32  |   0   |   26  |   0          |   -0.635                    |
| 31  |   0   |   32  |   0          |   -0.635                    |
| 36  |   1   |   31  |   1          |   -0.635                    |
| …   |   …   |   …   |   …          |   …                         |

 

 

The phenotype file provides information on the disease and covariates of
the study subjects in a tabular (row-column) format. Each row contains
space or tab delimited data specific to an individual. Variable names
should be specified in the first line of the file. The disease variable
should be listed first and can be followed by an arbitrary number of
covariates (or no covariate). In a case-control study, the disease
variable should be coded 0/1 to represent unaffected/affected. In a
cohort study which has two disease variables, the time variable should
be listed first and the indicator of disease second. Missing disease
variables or covariates are denoted as ‘.’.

 

### **genotype file**

 

|            |            |     |       |                                   |
|------------|------------|-----|-------|-----------------------------------|
|            |            |     |       | Example of a genotype file        |
| rs16977020 |   54706569 |   0 |   A C |   1 2 2 2 2 1 2 2 0 2 1 1 1 2 2 … |
| rs12903336 |   54715530 |   0 |   A G |   0 1 2 2 1 1 1 1 0 2 1 1 1 0 2 … |
| rs28678122 |   54743606 |   0 |   A C |   1 1 2 2 1 2 1 1 2 2 2 2 2 0 2 … |
| …          |   …        |   … |   …   |   …                               |

 

 

The genotype file provides genotype information for the study subjects
in a tabular (row-column) format. Each row contains space or tab
delimited data specific to a SNP. The columns follow the format

>   
> SNP_id position strand_orientation nucleo1 nucleo2 geno_1 … geno_n  

If the strand orientation information is not available, all
strand_orientation fields should be shown as 0. If this information is
available, flag 1 in the field indicates that the strand orientation in
the study data is different from the external panel (so the allele
coding of the external panel will be switched by the program) and flag 0
indicates strand consistency. In particular, if all the genotypes in the
external panel are in forward strand, then flag 1 means that the SNP in
the study was recorded on reverse strand. The strand orientation
information is only required for C/G and A/T SNPs. For all the other
types of SNPs, this field can be left 0. In nucleo1 and nucleo2 fields
are the nucleotides of the SNP and should be in the alphabetical order.
The genotypes are coded with 0, 1 and 2, referring to the count of
nucleo1. Missing genotype should be coded as 9.

 

### **external genotype file**

 

|            |            |       |                                     |
|------------|------------|-------|-------------------------------------|
|            |            |       | Example of an external file         |
| rs4774891  |   54807077 |   C T |   1 1 2 2 1 2 2 2 2 1 2 2 1 2 2 2 … |
| rs8025391  |   54808154 |   A T |   1 1 2 0 1 1 2 1 1 0 2 1 1 1 1 2 … |
| rs10518872 |   54809475 |   G T |   2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 … |
| …          |   …        |   …   |   …                                 |

 

 

The external genotype file follows the same format as the genotype file
except that the strand_orientation column is absent. position should be
in the same ascending or descending order as in the study genotype file.
Trio data should be entered first, followed by duos and unrelated
individuals. Trios and duos should be entered in family blocks. Within
each trio, the child genotype is entered last.

 

### **haplotype file**

 

<table style="width:89%;">
<colgroup>
<col style="width: 18%" />
<col style="width: 70%" />
</colgroup>
<tbody>
<tr class="odd">
<td></td>
<td>Example of a phased haplotype file</td>
</tr>
<tr class="even">
<td>rs4774891</td>
<td>  1 1 1 0 0 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1<br />
…</td>
</tr>
<tr class="odd">
<td>rs8025391</td>
<td>  0 0 1 0 0 1 1 0 1 1 1 1 1 1 0 1 0 1 0 1 1 1 1<br />
…</td>
</tr>
<tr class="even">
<td>rs10518872</td>
<td>  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1<br />
…</td>
</tr>
<tr class="odd">
<td>…</td>
<td>  …</td>
</tr>
</tbody>
</table>

 

 

The haplotype file for the external panel can be incorporated by using
the flag **-h**. The file format is

>   
> SNP_id phase1_1 phase2_1 phase1_2 phase2_2 …phase1_n phase2_n  

SNP_id should be in the same order as in the external genotype file.
Each subject contributes two columns (phase1_i phase2_i, i=1, …, n) with
0/1 coding, referring to the count of nucleo1 as in the external
genotype file. Typically, only phasing information on founders (mothers
and fathers) is provided for trios.

## **OUTPUT**

The file format is

> Typed SNP-id Position 1/0 MD Freq Estimate StdErr Z-Stat p-Value

For each untyped SNP, Typed is set to be “no”, 1/0 indicates the
nucleotides coded as 1 and 0, MD is the MD measure between the SNP and
the set of typed SNPs with the best prediction, Freq is the frequency of
1-coded allele in the external panel. Any untyped SNP with allele
frequency 0.0 or 1.0 in the external panel is excluded from the
analysis.

For each genotyped SNP, Typed is the proportion of non-missing
genotypes, MD is set to be ’-’, Freq is the frequency of 1-coded allele
in the study. Any SNP with allele frequency 0.0 or 1.0 in the study is
excluded from the analysis. The results for alleles with very low
minor-allele frequencies may not be stable and should be viewed with
great caution, especially for untyped SNPs or typed SNPs with
substantial missingness.

## **EXAMPLE**

The example includes a specification file “`GAWspec.txt`“, a phenotype
file “`GAWpheno.dat`“, a genotype file “`GAWgeno.dat`“, an external file
“`GAWexter.dat`“, and a haplotype file “`GAWhaplo.dat`“.

Enter the command

> \$ SNPMStat -sfile GAWspec.txt -pfile GAWpheno.dat -gfile GAWgeno.dat
> -efile GAWexter.dat -h -hfile GAWhaplo.dat -speed -ofile GAW

to obtain the results given in “`GAW_G.out`” and “`GAW_CPD$G.out`“.

## **REFERENCE**

Hu, Y. J., Lin, D. Y. and Zeng, D. (2010), “A General Framework for
Studying Genetic Effects and Gene-Environment Interactions with Missing
Data”, *Biostatistics*, in press.

Lin, D. Y., Hu, Y. and Huang, B. E. (2008), “Simple and Efficient
Analysis of SNP-Disease Association with Missing Genotype Data”,
*American Journal of Human Genetics*, 82, 444-452.

## **DOWNLOAD**

#### **SNPMStat for Linux \[updated July 13 2010\]**

executable (zip archive) **»**
[SNPMStat-4.0-linux.zip](http://dlin.web.unc.edu/wp-content/uploads/sites/1568/2013/01/SNPMStat-4.0-linux.zip)

#### **Example files \[updated July 13 2010\]**

zip archive **»**
[SNPMStat-4.0-example.zip](http://dlin.web.unc.edu/wp-content/uploads/sites/1568/2013/01/SNPMStat-4.0-example.zip)

## **VERSION HISTORY**

 

<table style="width:100%;">
<colgroup>
<col style="width: 4%" />
<col style="width: 7%" />
<col style="width: 86%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Version</th>
<th style="text-align: left;">Date</th>
<th style="text-align: left;">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">1.0</td>
<td style="text-align: left;">Oct. 2007</td>
<td style="text-align: left;">First version released</td>
</tr>
<tr class="even">
<td style="text-align: left;">1.1</td>
<td style="text-align: left;">May 8, 2008</td>
<td style="text-align: left;"><p><strong>Bug Fix:</strong></p>
<ul>
<li>Fixed bug in the logistic regression analysis. This bug only
affected the results of typed SNPs without missing values.</li>
</ul></td>
</tr>
<tr class="odd">
<td style="text-align: left;">2.0</td>
<td style="text-align: left;">Jul. 9, 2008</td>
<td style="text-align: left;"><p><strong>New Features:</strong></p>
<ul>
<li><p>1) Requires additional columns in input data files to include
nucleotide information and strand information.</p></li>
<li><p>2) Added internal check of the strand orientation between
reference panel and case-control data files. Added
<strong>-no_fix</strong> and <strong>-no_remove</strong> to control the
internal check.</p></li>
<li><p>3) Revised the example.</p></li>
<li><p>4) Added supplementary program SNPMStat_HM to convert the
reference panel from HapMap database to the required format.</p></li>
</ul></td>
</tr>
<tr class="even">
<td style="text-align: left;">2.1</td>
<td style="text-align: left;">Sep. 29, 2008</td>
<td style="text-align: left;"><p><strong>Bug Fix:</strong></p>
<ul>
<li>Fixed bug in reading <em>case_control.dat</em> when there is no
reference panel available. This bug only affected the results when
option <strong>-nr</strong> was specified.</li>
</ul></td>
</tr>
<tr class="odd">
<td style="text-align: left;">3.0</td>
<td style="text-align: left;">Oct. 14, 2008</td>
<td style="text-align: left;"><ul>
<li><p><strong>Bug fix:</strong></p>
<ul>
<li>Included SNPs in <em>case_control.dat</em> but not in
<em>reference.dat</em> into the analysis.</li>
</ul></li>
<li><p><strong>New Features:</strong></p>
<ul>
<li><p>1) Added supplementary program SNPMStat_CC to convert the
case-control data from PLINK format to the required format.</p></li>
<li><p>2) Required the physical position column for
<em>case_control.dat</em></p></li>
<li><p>3) Changed the format of <em>case_control.dat</em> when there is
no reference panel to be the same as the one with reference
panel.</p></li>
<li><p>4) Relaxed the requirement of rs number to allow any SNP
identifier.</p></li>
<li><p>5) Added a column with nucleotide coding information in the
output file.</p></li>
<li><p>6) Added the option <strong>-dom</strong> to allow the analysis
of dominant effect.</p></li>
<li><p>7) Added the option <strong>-speed</strong>.</p></li>
</ul></li>
</ul></td>
</tr>
<tr class="even">
<td style="text-align: left;">3.1</td>
<td style="text-align: left;">Dec. 17, 2008</td>
<td style="text-align: left;"><ul>
<li><p><strong>New Features:</strong></p>
<ul>
<li><p>1) Added the option <strong>-impute</strong> to allow the
imputation of untyped SNPs or missing values of typed SNPs.</p></li>
<li><p>2) Added the option <strong>-out_imp</strong> to specify the
output file for imputed genotypes.</p></li>
<li><p>3) Added the option <strong>-notest</strong> to allow the
suppression of association analysis.</p></li>
</ul></li>
</ul></td>
</tr>
<tr class="odd">
<td style="text-align: left;">3.2</td>
<td style="text-align: left;">Feb. 17, 2009</td>
<td style="text-align: left;"><ul>
<li><p><strong>Bug Fix:</strong></p>
<ul>
<li>Fixed a bug that may crash the program in certain cases.</li>
</ul></li>
<li><p><strong>New Feature:</strong></p>
<ul>
<li>with “-speed”, typed SNPs that are not in <em>phase.dat</em> are
analyzed by complete-case Armitage test and are not imputed even when
“-impute” is specified.</li>
</ul></li>
</ul></td>
</tr>
<tr class="even">
<td style="text-align: left;">4.0</td>
<td style="text-align: left;">Jul. 13, 2010</td>
<td style="text-align: left;"><ul>
<li><p><strong>New Feature:</strong></p>
<ul>
<li><p>1) Expanded the program to allow environmental factors and
gene-environment interactions. The environmental factors are allowed to
be correlated with genetic factors.</p></li>
<li><p>2) Expanded the program to allow cross-sectional and cohort
studies.</p></li>
<li><p>3) New format of input data fil</p></li>
</ul></li>
</ul></td>
</tr>
</tbody>
</table>

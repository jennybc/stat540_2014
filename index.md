---
css.file: _fileherd/markdown7.css
dirs.to.publish:
  - .
---

# STAT540 Statistical methods for high dimensional biology

## Course Information

### Credits and cross-listing

STAT 540 is a 3 credit course with a mandatory computing seminar.

Cross-listed as

  * STAT 540
  * BIOF 540
  * GSAT 540

### Instructors

* [Jennifer Bryan](http://www.stat.ubc.ca/~jenny), Course coordinator
  (email is jenny@stat.ubc.ca), Statistics and Michael Smith Labs

* [Gabriela Cohen-Freue](http://www.stat.ubc.ca/~gcohen/), Statistics

* [Paul Pavlidis](http://www.chibi.ubc.ca/faculty/pavlidis), CHiBi and Psychiatry

### TA(s)

* [Shaun Jackman](http://sjackman.ca), sjackman@gmail.com

* [Luolan (Gloria) Li](http://www.msl.ubc.ca/users/gloriali), glorialuolanli@gmail.com

### Google Group for Q & A (TAs will add students to group in due course)

[STAT540_2014](https://groups.google.com/forum/#!forum/stat540_2014)

### Time and Location

06 January 2014 - 07 April 2014

#### Lecture (Sec 201)
Time : Mon Wed 9:30 - 11am

Location : [ESB 4192](http://www.students.ubc.ca/classroomservices/buildings-and-classrooms/?code=ESB&room=4192)

#### Seminar / computing lab (S2A) -- REGISTRATION IS REQUIRED!
Time: *officially* runs Wed 12pm - 1pm; *unofficially* students are welcome to come after class around 11am and begin a ~1 hour guided analysis with TA support; TA will remain in the lab until 1pm to help those who start as 12pm and for general office hours.

Location: ESB 1042 __and 1046__

### Prerequisites.

Officially none BUT here in reality ...

* Statistics: you should have already taken university level introductory statistics course.

* Biology: No requirements, but you are expected to learn things like the difference between a DNA and RNA and a gene and a genome.

* R: no experience required but be prepared to do *a lot* of self-guided learning. Go ahead and start now by [installing R](http://www.r-project.org/) and the HIGHLY RECOMMENDED "integrated development environment" (IDE) [RStudio](http://www.rstudio.com)! Students are expected to run R on their own computer or a computer they have plenty of access to and control over. The best set-up, if possible, is to bring your own laptop to the computing seminars.

### Evaluation

* Homework. Two assignments worth 25 points each. Homework #1 will be assigned in early February and will be due in late February. Homework #2 will be assigned in early March and will be due in late March. *New for 2014: we may suggest/require that incremental progress be submitted along the way.*

* [Group project](project/index.html). Groups formed and projects conceived during January/February. Primary deliverable is a poster, presented in last class meeting. Each student also produces a short report. 40 points. More more information, [go here](project/index.html).

* 10 points for "other", e.g. participation in class, seminars, and the Google group, engagement with *small* computing exercises.

## Syllabus

### Week 1

seminar 00 | [R, RStudio Set Up & Basics](http://www.stat.ubc.ca/~jenny/STAT545A/block00_setup.html), borrowed from STAT 545A. Students complete/attempt on their own in advance. Bring any difficulties to first seminar.

lecture 01 | Introduction to high dimensional biology and the course (PP) | Mon Jan 06 | [slides as PDF](pvt/lect01_course-intro.pdf)

lecture 02 | Overview / review of probability and statistical inference, 1 of 2 (JB) | Wed Jan 08 | [slides as PDF](pvt/lect02_introToStatInf-probBasics.pdf)

seminar 01 | [R basics and exploring a small gene expression dataset](#seminar01) | Wed Jan 08

  * [R stuff](#seminar01) 11am - 12pm (or later, if necessary)
  * Molecular biology/genetics 101 (LL), 12pm - 1pm | [slides as PDF](pvt/lect00_biology-intro.pdf)

### Week 2

lecture 03 | Overview / review of probability and statistical inference, 2 of 2 (JB) | Mon Jan 13

lecture 04 | Exploratory analysis (PP) | Wed Jan 15

seminar 02 | Introduction To Simulation (LL) | Wed Jan 15

### Week 3

lecture 05 | Data QC and preprocessing (GC-F) | Mon Jan 20

lecture 06 | Statistical inference: two group comparisons, e.g. differential expression analysis (JB) | Wed Jan 22

seminar 03 | Introduction to R graphics (SJ) | Wed Jan 22

Fri Jan 24: Project groups should be formed. 

### Week 4

lecture 07 | Statistical inference: more than two groups --> linear models (JB prep, GCF deliver) | Mon Jan 27

lecture 08 | Statistical inference: linear models with 2 categorical covariates, greatest hits of linear models inference (JB prep, GCF deliver) | Wed Jan 29

seminar 04 | Data aggregation and two group testing (LL) | Wed Jan 29

Fri Jan 31: Initial project proposals due.

### Week 5

lecture 09 | Statistical inference: linear models including a quantitative covariate, fitting many linear models at once (JB) | Mon Feb 03

lecture 10 | Large scale inference: Empirical Bayes, limma (JB) | Wed Feb 05

seminar 05 | Fitting and interpretting linear models (low volume) (SJ) | Wed Feb 05

Fri Feb 07: Homework #1 will be assigned.

### Week 6

Mon Feb 10 is Family Day; no class

lecture 11 | Large scale inference: multiple testing (JB) | Wed Feb 12

seminar 06 | Fitting and interpretting linear models (high volume) (SJ), limma package (JB) | Wed Feb 12

Fri Feb 14: feedback to groups re: initial project proposals. Each group will be assigned an instructor or TA + instructor pair for extra support.

### Week 7

(Mon Feb 17 UBC closed for mid-term break)

(Wed Feb 19 UBC closed for mid-term break)

### Week 8

lecture 12 | Analysis of RNA-Seq data (PP), 1 of 2 | Mon Feb 24

Tues Feb 25: Homework #1 due 

lecture 13 | Analysis of RNA-Seq data (PP), 2 of 2 | Wed Feb 26

seminar 07 | RNA-Seq analysis (LL) | Wed Feb 26

Fri Feb 28: final project proposals due.

### Week 9

lecture 14 | Analysis of epigenetic data, focus on methylation (PP or designate) | Mon Mar 03
 
lecture 15 | Principal component analysis (PP) | Wed Mar 05

seminar 08 | Methylation analysis (LL) | Wed Mar 05

Fri Mar 07: Homework #2 assigned 

### Week 10

lecture 16 | Cluster analysis (GC-F) | Mon Mar 10

lecture 17 | Classification (GC-F) | Wed Mar 12

seminar 09 | Clustering and PCA (LL) | Wed Mar 12

### Week 11

lecture 18 | Model and variable selection: cross validation and regularization (GC-F) | Mon Mar 17

lecture 19 | Regularization cont'd, Proteomics and missingness (GC-F) | Wed Mar 19

seminar 10 | Supervised learning, cross validation, variable selection (SJ) | Wed Mar 19

### Week 12

lecture 20 | Enrichment analysis, gene networks, 1 of 2 (PP) | Mon Mar 24

lecture 21 | Enrichment analysis, gene networks, 2 of 2 (PP) | Wed Mar 26

seminar 11 | TA office hours during seminar time ... group project work | Wed Mar 26

Fi Mar 28: Homework #2 due.

### Week 13

lecture 22 | Resampling and the bootstrap (JB)| Mon Mar 31

lecture 23 | Guest lecture by STAT540 alum [Dr. Sohrab Shah](http://compbio.bccrc.ca/people/dr-sohrab-shah/) | Wed Apr 02

seminar 12 | TA office hours during seminar time ... group project work | Wed Apr 02

### Week 14

lecture 24 | Poster session 9:30am - 12:00pm | Wed Apr 09 __<-- NOTE THIS IS ON WEDNESDAY__, instead of Monday. Location: Room 101, aka the multi-purpose room, on ground floor of the [Michael Smith Building](http://www.maps.ubc.ca/PROD/index_detail.php?locat1=083).

## Seminars (guided analyses)

We will borrow some material from [STAT 545A Exploratory Analysis](http://www.stat.ubc.ca/~jenny/STAT545A/quick-index.html), in addition to using content specific to STAT 540.

  * seminar 00 | [R, RStudio Set Up & Basics](http://www.stat.ubc.ca/~jenny/STAT545A/block00_setup.html), borrowed from STAT 545A.
  * seminar 01 Wed Jan 08 <a id="seminar01"></a>
    - [Basics of R/RStudio, workspaces, and projects](http://www.stat.ubc.ca/~jenny/STAT545A/block01_basicsWorkspaceWorkingDirProject.html), borrowed from STAT 545A.
    - [Basic care and feeding of data in R](http://www.stat.ubc.ca/~jenny/STAT545A/block02_careFeedingData.html), borrowed from STAT 545A.
    - [R objects (beyond data.frames) and indexing](http://www.stat.ubc.ca/~jenny/STAT545A/block03_basicObjects.html), borrowed from STAT 545A.
    - [Explore a small gene expression dataset](seminars/seminar01_basic-data-exploration.html)
    - [Prep work for later use of Git, GitHub, Rpubs, etc.](seminars/seminar90_getting-started.html)
  * seminar 02 Wed Jan 15
    - [Probability and Simulation in R](seminars/seminar02_probability-simulation.html) <-- NOT HERE YET!
  * seminar 03 Wed Jan 22
    - Introduction to R graphics
  * seminar 04 Wed Jan 29
    - Data aggregation and two group testing
  * seminar 05 Wed Feb 05
    - Fitting and interpretting linear models (low volume)
  * seminar 06 Wed Feb 12
    - Fitting and interpretting linear models (high volume), limma package
  * seminar 07 Wed Feb 26
    - RNA-Seq data
  * seminar 08 Wed Mar 05
    - Methylation analysis
  * seminar 09 Wed Mar 12
    - Clustering and PCA
  * seminar 10 Wed Mar 19
    - Supervised learning, cross validation, variable selection

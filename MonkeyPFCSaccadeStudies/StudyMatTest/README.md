# Analysis of Saccades and Neuronal Activity

### Outline

1. Study Description
2. Setup
3. Preprocess
4. Analysis

### 1. Study Description

TODO

### 2. Setup

1. Follow the instructions common to all studies in the root [MonkeyStudies](https://github.com/SachsLab/MonkeyStudies) folder.

2. Indicate which sessions you want to analyze.
    - Use `MonkeyStudies/SessionInfo/SessionList.csv` as a guide.
    - Edit `MonkeyStudies/StudyLocationRule/Analysis/my_sessions.m` accordingly.

### 3. Preprocess

1. Change to `MonkeyStudies/StudyLocationRule/Anaysis/` as your working directory in Matlab.

2. Make sure the raw data are mounted.

3. Run `preproc_behav.m`. It will take a while.

4. Run `preproc_neural1.m`. This will convert the raw neural data to something suitable for mksort.

5. Do spike sorting.

6. Run `preproc_neural2.m`

You no longer need the raw data. You may now unmount it.

### 4. Analysis

1. Edit `MonkeyStudies/StudyLocationRule/Analysis/Setup/my_anaparams.m`.
	Here is where you define the analysis windows and the more general analysis parameters.
2. Run `MonkeyStudies/StudyLocationRule/Analysis/analyze_behav.m`.
	TODO:
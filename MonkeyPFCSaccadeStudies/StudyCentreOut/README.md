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
    - Use `./SessionInfo/SessionList.csv` as a guide.
    - Edit `./StudyCentreOut/Analysis/my_sessions.m` accordingly.

### 3. Preprocess

1. Change to `MonkeyStudies/StudyCentreOut/Anaysis/` as your working directory in Matlab.

2. Make sure the raw data are mounted.

3. Run `preproc_behav.m`. It will take a while.

4. Run `preproc_neural1.m`. This will convert the raw neural data to something suitable for mksort.

5. Do spike sorting with mksort.

6. Run `preproc_neural2.m`

You no longer need the raw data. You may now unmount it.

### 4. Analysis

1. Edit `./Setup/my_anaparams.m`.
	Here is where you define the analysis windows and the more general analysis parameters.
1. Run `./analyze_behav.m`.
	This produces figures showing the detected saccades, colour-coded according to the target location.
	You can also uncomment a section to plot individual trial behavioural data.
	This function does not produce output used in the paper, it is used only to verify the behaviour.
1. Run `./analyze_neural_epochs.m`.
	This calculates each neuron's tuning curve & PSTH, the class ANOVA for each task epoch, and the class ANOVA for the population for each task epoch.
	The result is saved to `../Output/Figures/tuning_epochs.mat`.
1. Run `./analyze_neural_timeseries.m`.
	This calculates the class ANOVA on firing rate for each neuron and for the population at each time-step in the trial. The result is saved to `../Output/mi_timeseries.mat`.
1. Run `./analyze_neural_classification.m`.
	This does a cross-validated classification analysis on multiple feature sets using multiple ML algorithms.
	This produces `../Output/classification.mat`
1. Run `./analyze_neural_ITR.m`.
	This analyzes the ITR as a function of trial-length using the neural trajectory and classifier with the best accuracy.
	This produces `../Output/ITR.mat`
1. Run `./analyze_neuron_adding.m`
	This analyzes the classification accuracy using 1 to N neurons using the best feature-set and ML algorithm.


### 5. Output

1. Run `./plot_exemplars.m`
	This generates candidate panels to be placed in Figure 3.
2. Run `./plot_neural_characterization.m`
	This creates `../Output/Figures/tuning_average.fig`, used for Figure 4, `mi_timeseries.fig` used for Figure 5, and `../Output/Table1.csv` used for Table1.
	It also prints some information to the command window that can be used in the manuscript.
1. Run `./plot_neural_char_traj.m`
	This uses a per-epoch population ANOVA to transform population activity into canonical trajectories, and compares their distance to GPFA trajectories.
	This creates Figures for each session. I chose session 10 (I think) for Figure 6. It also prints out some information.
4. Run `./plot_neural_classification.m`
	This creates `../Output/Figures/classification.fig`, used for Figure 7.
5. Run `./plot_neural_ITR.m`
	This creates `../Output/Figures/itr.fig` used for Figure 8.
	
### Notes on converting figures for publication

For Figures without image plots (i.e., only lines and symbols):

	Using Matlab, open the .fig and save as .svg
	Open the .svg in [Inkscape](https://inkscape.org/en/)
	I try to remove the background in Matlab, but check there is no background by double clicking on the image and trying to select and delete the background object.
	Change the ruler units to cm
	Select all items in the image, lock the H-V relative sizes, and change the horizontal size to match that of the journal. (e.g., for [JNS](http://www.jneurosci.org/site/misc/ifa_illustrations.xhtml#Figures); [JNP](http://www.the-aps.org/mm/Publications/Info-For-Authors/Preparing-Figures) )
	Click on File->Document Properties, under "Custom Size", expand Resize page to content and click on "Resize page to drawing or selection".
	Save the svg

For Figures with images

	Using Matlab, open the .fig and save as .tiff
	Open the .tiff in [GIMP](http://www.gimp.org/)
	Use the Image->Scale Image.. menus to set the size and resolution of the image.
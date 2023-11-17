# SDD_dots readme v0.8.1
Barcode and dot extractor for DNA stretched on glass

# How to use
With GUI:
Run `sdd_gui`. The GUI asks for a number of parameters and then runs `sdd_process_folder` routine which allows for analysis of large data sets.
Information about some specific settings can be found below.
--
Without GUI:
```
datafold = 'testfolder/';
[output,hPanelResult,images,movies,barcodes]  = sdd_script('sdd_settings.txt',[],datafold});
```
The following settings are available
```
130     | Pixel size(nm)                        | pxnm
300     | Width of LoG filter (nm)              | logSigmaNm
C=0     | Molecule image flag                   | barFlag
C=1     | Dots image flag                       | dotFlag
0       | Minimum log(EdgeScore)                | lowLim
10      | Minimum DotScore                      | dotScoreMin
1       | Minimum width (px)                    | widthLims(1)
Inf     | Maximum width (px)                    | widthLims(2)
50      | Minimum length (px)                   | lengthLims(1)
Inf     | Maximum length (px)                   | lengthLims(2)
2       | Edge margin for dots                  | dotMargin
0.8     | Minimum molecule eccentricity         | elim
0.4     | Minimum molecule-to-convex-hull ratio | ratlim
1       | Show score histograms                 | showScores
1       | Show detected molecules               | showMolecules
0       | Save detected molecules               | saveMolecules
0       | Save barcodes and dots                | saveBars
1       | Auto-threshold EdgeScore              | autoThreshBars
1       | Auto-threshold DotScore               | autoThreshDots
1       | Spline dot-detection                  | extractionMethod
5       | numSigmasAutoThresh                   | numSigmasAutoThresh
randomBars        | autoThreshDotsMethod (options - randomBars,meanstd) | autoThreshDotsMethod
0       | Remove non-uniform noise from dot images | denoiseDotImages
20     | length random barcode (for autothresh) | lenRandBar
0 |Set lower limit for number of standard deviations in molecule intensity from the background | sigmaBgLim
3 | Minimum distance (in pixels) from image edge for a molecule to be included in the analysis.  | edgeMargin
1 | Number of sigma_psf uncertainty for extract_barcodes.| deltaCut
0 | | showDotPeaks
```


As this software is set up to analyse pairs of images, one with DNA molecules and one with "dots", it requires particular naming of the images in the folder:
- The names of the images showing barcodes due to fluouresence by YOYO must contain the flag defined in `experiment.barFlag`.
- The names of the images showing dots must contain the flag defined in `experiment.dotFlag`.
Example: If the YOYO fluoresence image always contains "CH1", such that all YOYO images are on the form "XXCH1YY.tif", set `experiment.barFlag = 'CH1'`.
Similarly, if the dot images are named on the form "XXCH2YY.tif", set `experiment.dotFlag = 'CH2'`.

Update(October 2020): It is now possible to extract both barcodes and dots from single images, although this is strongly discouraged.
Set `experiment.barFlag` and `experiment.dotFlag` to empty character arrays to achieve this functionality.


One may tune the parameters which filter bad molecules, by changing the "User default settings" in the 'sdd_settings.txt' file or through SDD_Gui.
If `showScores=1`, the routine automatically plots a histogram of the scores generated by the regions in the images as a guide.
Additionally the detected molecules are numbered in the original image and the LoG filtered image are shown.

An option for automated choice of the thresholds for molecule and dot scores has been added.
To enforce automated thresholding of these, set `autoThreshBars=1` and/or `autoThreshDots=1`.

If `saveMolecules=1`, closeups of the detected molecules along with their barcode, if a such can be extracted, is saved in a "molecules" folder inside the target folder.
These should be readily insertable into the `CBC_Gui` routine. Similarly, if `saveBars=1`, the molecule barcodes and barcode from the dot image are saved in the folders "barcodes" and "dotbars" respectively.

For each image where dots have been detected, the number of barcodes, total APPARENT length of barcodes, number of detected dots(within the barcode regions) and the average dots per micron are printed in a `...results(I).txt` file, where (I) is an integer.
The program automatically increases (I) by 1 for each analysis if an old `...results(I).txt` file is present in the folder. When 10 such files are present, the new results are always saved as `...results10.txt`.

The algorithm scans through all levels of subfolder in the given directory and runs the `dnarec_skel` routine on all SUITABLE subfolders;
A subfolder is deemed "suitable" if
- Its name does not contain the phrase "molsandbars" or "movies"
- It has no subfolders inside it that fullfill the first point.
This means that if a folder has data images but also a subfolder which is does not contain the phrase "molsandbars" or "movies", the folder will be ignored.

Please cite this software as

[![DOI](https://zenodo.org/badge/547354110.svg)](https://zenodo.org/doi/10.5281/zenodo.10149648)




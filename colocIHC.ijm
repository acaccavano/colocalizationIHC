setBatchMode(true);
roiManager("deselect");
run("Select None");

closeTempImgs = true;  // Usually want this true, only disable for debug
verboseOption = false; // Displays (many) messages to the log
subBGOption   = true;  // Performs background subtractionb
subBGSize     = 50;    // Background substraction pixel size (Default = 20)

// Determine image names and split channels
rawName  = getTitle();
Stack.getDimensions(w, h, nChannel, nSlice, nFrame);
if (verboseOption) print("Width: " + w + ", Height: " + h + ", nChannel: " + nChannel + ", nSlice: " + nSlice + ", nFrame: " + nFrame);
Stack.getActiveChannels(actChStr);
chName   = newArray(nChannel);
tempImg  = newArray(nChannel);
actCh    = newArray(nChannel);
run("Duplicate...", "duplicate");
dupImg   = getTitle();
run("Split Channels");

// Determine channel names, active status, and background subtraction:
for (ch = 0; ch < nChannel; ch++) {
	chName[ch] = "C" + (ch+1) + "-" + dupImg;
	actCh[ch] = parseFloat(actChStr.substring(ch, ch+1));
	
	// Background subtraction:
	if (subBGOption && actCh[ch]) {
		selectWindow(chName[ch]);
		if (nSlice > 1) run("Subtract Background...", "rolling=" + subBGSize + " stack");
		else run("Subtract Background...", "rolling=" + subBGSize);
	}
}

// Determine number and type of ROIs
nROI = RoiManager.size;
if (verboseOption) print("nROI: " + nROI);

typeROI     = newArray(nROI);
nameROI     = newArray(nROI);
regionIdx   = newArray(nROI);
cellTypeIdx = newArray(nROI);
cellTypeCh  = newArray(nROI);

nRegion   = 0;
nCellType = 0;

// Determine ROI type:
for (i = 0; i < nROI; i++) {
	RoiManager.select(i);
	typeROI[i] = Roi.getType;
	nameROI[i] = Roi.getName;
	if (typeROI[i] == 'polygon') {
		regionIdx[nRegion] = i;
		if (verboseOption) print("regionIdx[" + nRegion + "] = " + regionIdx[nRegion] + " (" + nameROI[i] + ")");
		nRegion++;
	} else if (typeROI[i] == 'point') {
		cellTypeIdx[nCellType] = i;
		// Special logic to assign channel of cell
		if (nameROI[i].indexOf("-") == 1) {
			cellTypeCh[nCellType] = parseFloat(nameROI[i].substring(0, nameROI[i].indexOf("-")));
			nameROI[i] = nameROI[i].substring(nameROI[i].indexOf("-") + 1);
			if (verboseOption) print("cellTypeCh[" + nCellType + "]: " + cellTypeCh[nCellType] + " (" + nameROI[i] + ")");
		}
		if (verboseOption) print("cellTypeIdx[" + nCellType + "]: " + cellTypeIdx[nCellType] + " (" + nameROI[i] + ")");
		nCellType++;
	}
}

regionIdx   = Array.trim(regionIdx, nRegion);
cellTypeIdx = Array.trim(cellTypeIdx, nCellType);
cellTypeCh  = Array.trim(cellTypeCh, nCellType);
nCell       = newArray(nCellType);

outSlice  = newArray(nRegion * nSlice);
outRegion = newArray(nRegion * nSlice);
outName   = newArray(nRegion * nSlice);
outArea   = newArray(nRegion * nSlice);

// Initialize mean and stdev variables (max of 4 channels possible)
for (ch = 0; ch < nChannel; ch++) {
	if (ch==0) {
		outMnCh1 = newArray(nRegion * nSlice);
		outSDCh1 = newArray(nRegion * nSlice);
	} else if (ch==1) {
		outMnCh2 = newArray(nRegion * nSlice);
		outSDCh2 = newArray(nRegion * nSlice);
	} else if (ch==2) {
		outMnCh3 = newArray(nRegion * nSlice);
		outSDCh3 = newArray(nRegion * nSlice);
	} else if (ch==3) {
		outMnCh4 = newArray(nRegion * nSlice);
		outSDCh4 = newArray(nRegion * nSlice);
	}
}

// Determine number of unique channel combinations:
nCorr = 0;
for (ch1 = 0; ch1 < nChannel; ch1++) {
	for (ch2 = ch1; ch2 < nChannel; ch2++) {
		if (actCh[ch1] && actCh[ch2]) {
			if (ch1!=ch2) nCorr++;
		}
	}
}

corrCh1  = newArray(nCorr);
corrCh2  = newArray(nCorr);
corrName = newArray(nCorr);

// Assign variables for unique channel combinations (max possible of 6 combos for 4 channel images)
nCorr = 0;
for (ch1 = 0; ch1 < nChannel; ch1++) {
	for (ch2 = ch1; ch2 < nChannel; ch2++) {
		if (actCh[ch1] && actCh[ch2]) {
			if (ch1!=ch2) {
				corrCh1[nCorr] = ch1;
				corrCh2[nCorr] = ch2;
				corrName[nCorr] = "PCorr_Ch" + (ch1+1) + "-Ch" + (ch2+1);
				if (nCorr==0)      pCorr1 = newArray(nRegion * nSlice);
				else if (nCorr==1) pCorr2 = newArray(nRegion * nSlice);
				else if (nCorr==2) pCorr3 = newArray(nRegion * nSlice);
				else if (nCorr==3) pCorr4 = newArray(nRegion * nSlice);
				else if (nCorr==4) pCorr5 = newArray(nRegion * nSlice);
				else if (nCorr==5) pCorr6 = newArray(nRegion * nSlice);
				nCorr++;
			}
		}
	}
}

// Determine Region & Area
idx = 0;
for (i = 0; i < nRegion; i++) {
	
	// Calculate region area:
	RoiManager.select(regionIdx[i]);
	run("Set Measurements...", "area redirect=None decimal=6");
	roiManager("multi-measure measure_all");
	
	for (sl = 0; sl < nSlice; sl++) {

		// Update region names and area:
		outSlice[idx]  = sl + 1;
		outRegion[idx] = i + 1;
		outName[idx]   = nameROI[regionIdx[i]];
		outArea[idx]   = getResult("Area", sl)*0.000001;
		if (verboseOption) print(idx + ": Area (" + outName[idx] +  ") = " + outArea[idx] + " mm^2");
		idx++;
	}
	run("Clear Results");
}

// Determine Mean, SD, & Pearson's Correlation for each region and channel combination:
for (i = 0; i < nRegion; i++) {
	for (ch = 0; ch < nChannel; ch++) {
		selectWindow(chName[ch]);

		// Calculate Mean and SD for active channel and slice
		if (actCh[ch]) {
			
			RoiManager.select(regionIdx[i]);
			run("Set Measurements...", "mean standard redirect=None decimal=6");
			roiManager("multi-measure measure_all");
			
			idx = i * nSlice;
			for (sl = 0; sl < nSlice; sl++) {
				// Assign to array (works for up to 3 channels max):
				if (ch==0) {
					outMnCh1[idx] = getResult("Mean", sl);
					outSDCh1[idx] = getResult("StdDev", sl);
					if (verboseOption) print(idx + ": Mean_Ch" + (ch+1) + " (" + outName[idx] +  ") = " + outMnCh1[idx] + " (" + outSDCh1[idx] + " SD)");
				} else if (ch==1) {
					outMnCh2[idx] = getResult("Mean", sl);
					outSDCh2[idx] = getResult("StdDev", sl);
					if (verboseOption) print(idx + ": Mean_Ch" + (ch+1) + " (" + outName[idx] +  ") = " + outMnCh2[idx] + " (" + outSDCh2[idx] + " SD)");
				} else if (ch==2) {
					outMnCh3[idx] = getResult("Mean", sl);
					outSDCh3[idx] = getResult("StdDev", sl);
					if (verboseOption) print(idx + ": Mean_Ch" + (ch+1) + " (" + outName[idx] +  ") = " + outMnCh3[idx] + " (" + outSDCh3[idx] + " SD)");
				} else if (ch==3) {
					outMnCh4[idx] = getResult("Mean", sl);
					outSDCh4[idx] = getResult("StdDev", sl);
					if (verboseOption) print(idx + ": Mean_Ch" + (ch+1) + " (" + outName[idx] +  ") = " + outMnCh4[idx] + " (" + outSDCh4[idx] + " SD)");
				}
				idx++;
			}
			run("Clear Results");	
		}
	
		// Create new duplicate image normalized by mean and SD
		// Must create for all channels (even non-active) to avoid subsequent bug
		if (nSlice > 1) run("Duplicate...", "ignore duplicate");
		else run("Duplicate...", "ignore");
		tempImg[ch] = getTitle();
		
		if (actCh[ch]) {
			run("32-bit");
			idx = i * nSlice;
			for (sl = 0; sl < nSlice; sl++) {
				RoiManager.select(regionIdx[i]); // Selecting ROI can change slice, so make sure to update
				Stack.setSlice(sl+1);
				
				subStr = "value=";
				divStr = "value=";
				
				if (ch==0) {
					subStr = subStr + outMnCh1[idx];
					divStr = divStr + outSDCh1[idx];
				} else if (ch==1) {
					subStr = subStr + outMnCh2[idx];
					divStr = divStr + outSDCh2[idx];
				} else if (ch==2) {
					subStr = subStr + outMnCh3[idx];
					divStr = divStr + outSDCh3[idx];
				} else if (ch==3) {
					subStr = subStr + outMnCh4[idx];
					divStr = divStr + outSDCh4[idx];
				}
				
				if (nSlice > 1) {
					subStr = subStr + " slice";
					divStr = divStr + " slice";
				}
				
				run("Subtract...", subStr);
				run("Divide...", divStr);
				
				// Crop outside ROI
				run("Make Inverse");
				run("Clear", "slice");
				resetMinAndMax();
				idx++;
			}
		}
	}

	// Calculate pair-wise Pearson's correlations
	for (j = 0; j < nCorr; j++) {
		if (nSlice > 1) imageCalculator("Multiply create 32-bit stack", tempImg[corrCh1[j]], tempImg[corrCh2[j]]);
		else imageCalculator("Multiply create 32-bit", tempImg[corrCh1[j]], tempImg[corrCh2[j]]);
		multImg = getTitle();
		RoiManager.select(regionIdx[i]);
		run("Set Measurements...", "mean redirect=None decimal=6");
		roiManager("multi-measure measure_all");

		idx = i * nSlice;
		for (sl = 0; sl < nSlice; sl++) {

			// Assign to array (works for up to 3 comparisons (3 channels max):
			if (j==0) {
				pCorr1[idx] = getResult("Mean", sl);
				if (verboseOption) print(idx + ": " + corrName[j] + " (" + outName[idx] +  ") = " + pCorr1[idx]);
			} else if (j==1) {
				pCorr2[idx] = getResult("Mean", sl);
				if (verboseOption) print(idx + ": " + corrName[j] + " (" + outName[idx] +  ") = " + pCorr2[idx]);
			} else if (j==2) {
				pCorr3[idx] = getResult("Mean", sl);
				if (verboseOption) print(idx + ": " + corrName[j] + " (" + outName[idx] +  ") = " + pCorr3[idx]);
			} else if (j==3) {
				pCorr4[idx] = getResult("Mean", sl);
				if (verboseOption) print(idx + ": " + corrName[j] + " (" + outName[idx] +  ") = " + pCorr4[idx]);
			} else if (j==4) {
				pCorr5[idx] = getResult("Mean", sl);
				if (verboseOption) print(idx + ": " + corrName[j] + " (" + outName[idx] +  ") = " + pCorr5[idx]);
			} else if (j==5) {
				pCorr6[idx] = getResult("Mean", sl);
				if (verboseOption) print(idx + ": " + corrName[j] + " (" + outName[idx] +  ") = " + pCorr6[idx]);
			}
			idx++;
		}
		run("Clear Results");

		// Close the multiplied Pearons's Corrleation Image:
		if (closeTempImgs) {
			selectWindow(multImg);
			close();
		}
	}

	// Close all temp normalized images
	if (closeTempImgs) {
		for (ch = 0; ch < nChannel; ch++) {
			selectWindow(tempImg[ch]);
			close();
		}
	}
}

// Cell Counting and Mean Cell Intensities
for (j = 0; j < nCellType; j++) {
	
	// Initialize number arrays (works for up to 6 cell types):
	if (j==0)      nCellInc1 = newArray(nRegion);
	else if (j==1) nCellInc2 = newArray(nRegion);
	else if (j==2) nCellInc3 = newArray(nRegion);
	else if (j==3) nCellInc4 = newArray(nRegion);
	else if (j==4) nCellInc5 = newArray(nRegion);
	else if (j==5) nCellInc6 = newArray(nRegion);
	
	// Initialize mean arrays (works for up to 6 cell types) - only perform if cell assigned to channel
	if (cellTypeCh[j] > 0) {
		if (j==0)      mnCellInc1 = newArray(nRegion);
		else if (j==1) mnCellInc2 = newArray(nRegion);
		else if (j==2) mnCellInc3 = newArray(nRegion);
		else if (j==3) mnCellInc4 = newArray(nRegion);
		else if (j==4) mnCellInc5 = newArray(nRegion);
		else if (j==5) mnCellInc6 = newArray(nRegion);
	}
	
	// Determine number and position of cells:
	RoiManager.select(cellTypeIdx[j]);
	nCell[j] = Roi.size;
	Roi.getCoordinates(xPos, yPos);
	
	// Calculate mean intensity of each cell - only perform if cell assigned to channel
	if (cellTypeCh[j] > 0) {
		selectWindow(chName[cellTypeCh[j]-1]);
		run("Set Measurements...", "mean redirect=None decimal=6");
		roiManager("Measure");
		cellMean = newArray(nCell[j]);
		for (k = 0; k < nCell[j]; k++) cellMean[k] = getResult("Mean", k);
		run("Clear Results");
	}
	
	if (verboseOption) {
		print("nCell[" + j + "] (" + nameROI[cellTypeIdx[j]] + "): " + nCell[j]);
		if (cellTypeCh[j] > 0) {
			print("Cell:", "X:", "Y:", "Intensity:");
			for (k = 0; k < nCell[j]; k++) print(k, xPos[k], yPos[k], cellMean[k]);
		} else {
			print("Cell:", "X:", "Y:");
			for (k = 0; k < nCell[j]; k++) print(k, xPos[k], yPos[k]);
		}
	}
	
	// Loop over each region:
	for (i = 0; i < nRegion; i++) {
		RoiManager.select(regionIdx[i]);
		nCellInc = 0; // Number of cells included within the region
		cellMeanInc = newArray(nCell[j]);
		
		// Loop over each cell to obtain number of included cells
		for (k = 0; k < nCell[j]; k++) {
			if (Roi.contains(xPos[k], yPos[k])) {
				if (cellTypeCh[j] > 0) cellMeanInc[nCellInc] = cellMean[k];
				nCellInc++;
			}
		}

 		// Assign to correct arry:
		if (j==0)      nCellInc1[i] = nCellInc;
		else if (j==1) nCellInc2[i] = nCellInc;
		else if (j==2) nCellInc3[i] = nCellInc;
		else if (j==3) nCellInc4[i] = nCellInc;
		else if (j==4) nCellInc5[i] = nCellInc;
		else if (j==5) nCellInc6[i] = nCellInc;
		
		// If cell assigned to channel, determine average mean intensities:
		if (cellTypeCh[j] > 0) {

			cellMeanInc = Array.trim(cellMeanInc, nCellInc);
			Array.getStatistics(cellMeanInc, min, max, mnCellInc, std);
			
			// Assign to correct array:
			if (j==0)      mnCellInc1[i] = mnCellInc;
			else if (j==1) mnCellInc2[i] = mnCellInc;
			else if (j==2) mnCellInc3[i] = mnCellInc;
			else if (j==3) mnCellInc4[i] = mnCellInc;
			else if (j==4) mnCellInc5[i] = mnCellInc;
			else if (j==5) mnCellInc6[i] = mnCellInc;

		}

		if (verboseOption) {
			if (cellTypeCh[j] > 0) print("n" + nameROI[cellTypeIdx[j]] + " (" + outName[i] +  "): " + nCellInc + ", mean = " + mnCellInc);
			else print("n" + nameROI[cellTypeIdx[j]] + " (" + outName[i] +  "): " + nCellInc);
		}
	}
}

// Close duplicated images
if (closeTempImgs) {
	for (ch = 0; ch < nChannel; ch++) {
		selectWindow(chName[ch]);
		close();
	}
}

// Assign arrays to output table
outTable = substring(rawName, 0, lengthOf(rawName) - 4) + "_Output";
Table.create(outTable);
Table.setColumn("Region", outRegion);
Table.setColumn("Name", outName);
Table.setColumn("Slice", outSlice);
Table.setColumn("Area_mm2", outArea);

for (ch = 0; ch < nChannel; ch++) {
	if (actCh[ch]) {
		if (ch==0)      Table.setColumn("Mean_Ch" + (ch+1), outMnCh1);
		else if (ch==1) Table.setColumn("Mean_Ch" + (ch+1), outMnCh2);
		else if (ch==2) Table.setColumn("Mean_Ch" + (ch+1), outMnCh3);
		else if (ch==3) Table.setColumn("Mean_Ch" + (ch+1), outMnCh4);
	}
}

for (ch = 0; ch < nChannel; ch++) {
	if (actCh[ch]) {
		if (ch==0)      Table.setColumn("SD_Ch" + (ch+1), outSDCh1);
		else if (ch==1) Table.setColumn("SD_Ch" + (ch+1), outSDCh2);
		else if (ch==2) Table.setColumn("SD_Ch" + (ch+1), outSDCh3);
		else if (ch==3) Table.setColumn("SD_Ch" + (ch+1), outSDCh4);
	}
}

for (j = 0; j < nCorr; j++) {
	if (j==0)      Table.setColumn(corrName[j], pCorr1);
	else if (j==1) Table.setColumn(corrName[j], pCorr2);
	else if (j==2) Table.setColumn(corrName[j], pCorr3);
	else if (j==3) Table.setColumn(corrName[j], pCorr4);
	else if (j==4) Table.setColumn(corrName[j], pCorr5);
	else if (j==5) Table.setColumn(corrName[j], pCorr6);
}

for (j = 0; j < nCellType; j++) {
	if (j==0)      Table.setColumn("n" + nameROI[cellTypeIdx[j]], nCellInc1);
	else if (j==1) Table.setColumn("n" + nameROI[cellTypeIdx[j]], nCellInc2);
	else if (j==2) Table.setColumn("n" + nameROI[cellTypeIdx[j]], nCellInc3);
	else if (j==3) Table.setColumn("n" + nameROI[cellTypeIdx[j]], nCellInc4);
	else if (j==4) Table.setColumn("n" + nameROI[cellTypeIdx[j]], nCellInc5);
	else if (j==5) Table.setColumn("n" + nameROI[cellTypeIdx[j]], nCellInc6);
}

for (j = 0; j < nCellType; j++) {
	if (cellTypeCh[j] > 0) {
		if (j==0)      Table.setColumn("mean" + nameROI[cellTypeIdx[j]], mnCellInc1);
		else if (j==1) Table.setColumn("mean" + nameROI[cellTypeIdx[j]], mnCellInc2);
		else if (j==2) Table.setColumn("mean" + nameROI[cellTypeIdx[j]], mnCellInc3);
		else if (j==3) Table.setColumn("mean" + nameROI[cellTypeIdx[j]], mnCellInc4);
		else if (j==4) Table.setColumn("mean" + nameROI[cellTypeIdx[j]], mnCellInc5);
		else if (j==5) Table.setColumn("mean" + nameROI[cellTypeIdx[j]], mnCellInc6);
	}
}

// Save ouput table
saveAs("Results");

setBatchMode(false);
BatchMode = true;

GaussianBlurRad = 1.5;
NoiseTol = 8;

run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel global");
run("Set Measurements...", "  mean centroid median redirect=None decimal=2");
InitialStackID = getImageID();
InitialStackTitle = getTitle();	

// Copy stack + pre-filtering
run("Duplicate...", "title=Copy duplicate range=1-"+d2s(nSlices,0));
run("Gaussian Blur...", "sigma="+d2s(GaussianBlurRad,2)+" stack");
rename("PreProcessed");


// Batch mode
if(BatchMode==true)
{
	setBatchMode("exit & display");
	setBatchMode(true);
}

// Binary mask generation (apply seeded watershed to all slices)
newImage("ParticlesStack", "8-bit Black", getWidth(), getHeight(), nSlices);
for(i=1;i<=nSlices;i++)
{
	selectImage("PreProcessed");
	setSlice(i);
	run("Find Maxima...", "noise="+d2s(NoiseTol,2)+" output=[Segmented Particles] light");
	rename("Particles");
	run("Copy");
	selectImage("ParticlesStack");
	setSlice(i);
	run("Paste");
	selectImage("Particles");
	close();
}
selectImage("ParticlesStack");
run("Invert", "stack");
selectImage("PreProcessed");
close();
selectImage("ParticlesStack");
run("Select None");

run("Save", "save=/Users/eyal/Desktop/segmentedMovie.tif");

if(BatchMode==true)setBatchMode("exit & display");


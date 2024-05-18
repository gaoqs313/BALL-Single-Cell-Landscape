
SAMPLE=$1

## Step 1: normalize the input matrix and prepare the run parameters
cnmf prepare --output-dir cNMF --name $SAMPLE -c $SAMPLE.h5ad -k 3 4 5 6 7 8 9 10 --n-iter 100 --worker-index 0 --total-workers 1 --seed 14 --numgenes 2000

## Step 2: factorize the matrix
cnmf factorize --output-dir cNMF --name $SAMPLE --worker-index 0 --total-workers 1

## Step 3: combine the individual spectra results files for each K into a merged file
cnmf combine --output-dir cNMF --name $SAMPLE
rm cNMF/$SAMPLE/cnmf_tmp/$SAMPLE.spectra.k_*.iter_*.df.npz

## Step 4: select an optimal K by considering the trade-off between stability and error
cnmf k_selection_plot --output-dir cNMF --name $SAMPLE

## Step 5: set --local-density-threshold to be 2 for the first run after determine K
for K in {3..10}
do
	cnmf consensus --output-dir cNMF --name $SAMPLE --components $K --local-density-threshold 2 --show-clustering
done

## Step 6: set --local-density-threshold based on the histogram from previous step
for K in {3..10}
do
	for DT in 0.01 0.02 0.03 0.04 0.05
	do
		cnmf consensus --output-dir cNMF --name $SAMPLE --components $K --local-density-threshold $DT --show-clustering
	done
done


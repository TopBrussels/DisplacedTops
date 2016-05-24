array=( FullSamplesMuMuV0.xml FullSamplesElElV0.xml DisplacedTopsSignal.xml FullSamplesbbMuV0.xml FullSamplesbbElV0.xml )

for f in "${array[@]}"
do
    echo $f
    echo "about to run the xmlFiller!!!"
    python xmlFiller.py $f
done


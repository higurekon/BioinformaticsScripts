for cancer in {"CRC",}
do
mkdir $cancer
tar -xvf $cancer".tar" -C $cancer --wildcards 'FILE_SAMPLE_MAP.txt'
tar -xvf $cancer".tar" -C $cancer --wildcards '*nationwidechildrens.org_clinical_patient*'
tar -xvf $cancer".tar" -C $cancer --wildcards '*.isoforms.results'

prefix="allsamples_isoforms" 
output=$cancer/$prefix".txt"
filenames=$cancer/$prefix"_filenames.txt"

touch $output # Create output file
touch $filenames

awk '{print FILENAME;nextfile}' $cancer/{,**/}*isoforms.results > $filenames
awk '{a[$1]=a[$1] FS $3}END{for(i in a) print i,a[i]}' $cancer/{,**/}*isoforms.results > $output

sed -i 's/  /\t/g' $output
sed -i 's/ /\t/g' $output

done

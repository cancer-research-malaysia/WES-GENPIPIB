SLX_ID=$1; 
TUM_ID=$2; 

### trim
#ls |
#echo `date "+%Y-%m-%d %H:%M:%S"` "SLX-${SLX_ID} Trim Galore start" >> log.txt
file=$(echo $SLX_ID | cut -d '.' -f 1)
echo ${file} | xargs -P 4 -I %% sh -c "aws s3 ls crm.sequencing.raw.data.sharing/batch1/SLX-%%/" | grep SLX-${SLX_ID}. |
awk '{print $NF}'|
grep 'fq.gz$' | grep -v 0000 |
cut -d '.' -f 1,2,3,4 | 
sort | uniq | 
xargs -P 2 -I AAA sh -c  \
" echo AAA " ;
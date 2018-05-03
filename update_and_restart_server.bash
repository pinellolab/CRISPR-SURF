git pull 
docker stop $(docker ps -q --filter ancestor=pinellolab/crisprsurf)
docker run -p 9993:9993 \
-v /Users/sailor/Projects/CRISPR-SURF/SURF:/SURF \
-d -it pinellolab/crisprsurf SURF_webapp

#-v /Volumes/pinello/PROJECTS/2016_12_STREAM/STREAM_DB_top_15/:/STREAM/precomputed \

#-v /Volumes/Data/STREAM/UPLOADS_FOLDER:/tmp/UPLOADS_FOLDER \
#-v /Volumes/Data/STREAM/RESULTS_FOLDER:/tmp/RESULTS_FOLDER \

#-v /Volumes/Data/STREAM/tmp:/tmp \
#-v /Volumes/Data/STREAM_DB_top_15/:/STREAM/precomputed \


./prog NUM_PARTITIONS ENV_FILE_1 ENV_FILE_2 GEO_FILE_1 GEO_FILE_2 NUM_GPUS

# Example
Using 8 GPUs to partition sports and cemetery into 1024 parts

./prog 1024 ~/data/sports_env ~/data/cemetery_env ~/data/sports ~/data/cemetery 8


./prog 1024 ../../data/sports_data_env ../../data/cemet_data_env ../../data/sports_data ../../data/cemet_data 8

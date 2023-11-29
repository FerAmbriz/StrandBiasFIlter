#!/bin/bash

merge=$1
folder=$2

bash bin_size_run.sh $merge
bash run_bin.sh $folder

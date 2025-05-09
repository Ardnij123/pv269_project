#!/bin/bash

cd "$(dirname $0)"

conda install -n pv269_project $1 && \
    echo "conda install -n pv269_project $1" >> conda_create.sh

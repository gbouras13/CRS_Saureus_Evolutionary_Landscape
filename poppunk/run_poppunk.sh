#!/bin/bash

# assumes poppunk is installed in a conda env called poppunk

conda activate poppunk

poppunk_assign --db staphylococcus_aureus_v1_refs --query poppunk_input.txt \
--output poppunk_clusters --threads 8
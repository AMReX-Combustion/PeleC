#!/bin/bash -l

#Script that runs the nightly PeleC tests at NREL on Rhodes

set -e

cd /projects/ecp/combustion/pelec-testing2/logs && \
nice -n19 ionice -c3 \
/projects/ecp/combustion/pelec-testing2/pelec/Testing/test-pelec.sh &> \
"test-pelec-$(date +%Y-%m-%d).log"

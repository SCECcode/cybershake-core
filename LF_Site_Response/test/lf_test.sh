#!/bin/bash

ulimit -c unlimited

export MALLOC_CHECK_=2
../bin/lf_site_response seis_in=Seismogram_PAS_231_0_0_lf.grm seis_out=Seismogram_PAS_231_0_0_lf_site_response.grm slat=34.148426 slon=-118.17119 module=cb2014 vs30_model=cvmsi vs30=821.0

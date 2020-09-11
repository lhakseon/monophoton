#!/bin/sh
export pwd=/nfs_scratch/lhakseon/monophoton/CMSSW_10_2_3/src
cat>Job_${6}.sh<<EOF
#!/bin/sh
source /cvmfs/cms.cern.ch/cmsset_default.sh 
cd /nfs_scratch/lhakseon/monophoton/CMSSW_10_2_3/src
eval `scramv1 runtime -sh`
cd \${_CONDOR_SCRATCH_DIR}
./${1} ${2} ${3} ${4} ${5}
EOF

chmod 775 Job_${6}.sh

cat>condor_${6}<<EOF
x509userproxy = /tmp/x509up_u10053
executable = ./Job_${6}.sh
notification         = never
whenToTransferOutput = On_Exit
shouldTransferFiles  = yes
requirements = (TARGET.UidDomain == "hep.wisc.edu" && TARGET.HAS_CMS_HDFS)
on_exit_remove       = (ExitBySignal == FALSE && (ExitCode == 0 || ExitCode == 42 || NumJobStarts>3))
+IsFastQueueJob      = True
getenv = true
request_memory       = 1992
request_disk         = 2048000
transfer_input_files = ${1}
output               = \$(Cluster)_\$(Process)_${6}.out
error                = \$(Cluster)_\$(Process)_${6}.err
Log                  = \$(Cluster)_\$(Process)_${6}.log
Queue
EOF

condor_submit condor_${6}

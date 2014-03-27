# this script sets up the proper environment variables for a PBS job.

#PBS -V
#PBS -o stdout.log
#PBS -e stderr.log
##PBS -o /dev/null
##PBS -e /dev/null
#PBS -l mem=5gb

# go to directory from which job was submitted.  If this job is being run
# (not submitted to a queue) this should not run.
if [ ! -z $PBS_O_WORKDIR]; then
 echo "going to PBS working directory: ${PBS_O_WORKDIR}"
 cd ${PBS_O_WORKDIR}
fi

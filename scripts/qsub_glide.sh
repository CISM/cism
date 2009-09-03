# simple shell script to be used with SGE for launching glide
# tell SGE to change into current working directory
#$-cwd

echo Launching $1 on $HOSTNAME

$*

pushd  ~/projects/galaxy/galaxy/
./run.sh --stop-daemon
echo >paster.log
./run.sh --daemon
popd
tail -f  ~/projects/galaxy/galaxy/paster.log| grep -v "DEBUG"

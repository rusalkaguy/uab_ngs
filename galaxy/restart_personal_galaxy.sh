pushd  ~/projects/galaxy/galaxy/
./run.sh --stop-daemon
echo >paster.log
./run.sh --daemon
popd

echo "Waiting for galaxy to start..."
RC=1
while [ $RC == 1 ]; do
    echo "."
    sleep 1
    grep "^serving on http" paster.log 
    RC=$?
done
    

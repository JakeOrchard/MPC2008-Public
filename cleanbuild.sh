# shell script that deletes all the input and output folders in the subdirectories of the current directory

# find all subdirectories in the current directory
for dir in $(find . -maxdepth 1 -type d)
do
    # delete the input and output folders in the subdirectories
    rm -rf $dir/input $dir/output
done


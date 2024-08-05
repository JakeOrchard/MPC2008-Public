# takes as input a text file with a list of files to copy, creates a folder with the name of the text file if it does not yet exist, and copies the files listed in the text file to the folder with the name of the text file
# Usage: ./copytablesandfigures.sh <file>
# Example: ./copytablesandfigures.sh 2008tablesandfigures.txt

# check if the user has provided the correct number of arguments
if [ $# -ne 1 ]
then
    echo "Usage: ./copytablesandfigures.sh <file>"
    echo "Example: ./copytablesandfigures.sh 2008tablesandfigures.txt"
    exit 1
fi

# check if the file provided by the user exists
if [ ! -f $1 ]
then
    echo "The file $1 does not exist"
    exit 1
fi

# create a folder with the name of the file if it does not yet exist
foldername=../$(echo $1 | cut -d '.' -f 1)

if [ ! -d $foldername ]
then
    mkdir $foldername
fi

# if the folder exists, delete all the files in the folder
if [ "$(ls -A $foldername)" ]
then
    rm $foldername/*
fi

# copy the files listed in the text file to the folder with the name of the text file
while read line
do
    cp $line $foldername
done < $1

echo "Files copied successfully to $foldername"

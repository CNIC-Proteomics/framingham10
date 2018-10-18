INFILE="D:\projects\framingham_score\data\AWHS_params_framingham.xlsx"
OUFILE="D:\projects\framingham_score\data\AWHS_params_framingham.csv"

echo "** Calculate the Framingham 10"
echo "python calc_framingham10.py -i ${INFILE} -o ${OUFILE}"
python calc_framingham10.py -i ${INFILE} -o ${OUFILE}


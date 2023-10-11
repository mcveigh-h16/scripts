virtualenv myenv2 
source myenv2/bin/activate
pip install biopython
pip install pandas
pip install datetime

python processITS2.0.py protistITS2.gbk protistITS2.out

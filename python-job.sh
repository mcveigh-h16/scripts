virtualenv myenv 
source myenv/bin/activate.csh
pip install pandas
pip install biopython
python ParseCMscan1.65.py fullfungalnew.out

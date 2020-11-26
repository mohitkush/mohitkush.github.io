"""rdkit may not be avilable for some python users. If you are using Conda environment, there's no problem, it should work fine but if it doesn't here is the link to google colab where you can run rdkit easily """

import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
import requests
import json
import pickle


def mol_data(smiles):
    
    smiles['mol'] = smiles['CanonicalSMILES'].apply(lambda x: Chem.MolFromSmiles(x)) 
    smiles['mol_w'] = smiles['mol'].apply(lambda x: Descriptors.ExactMolWt(x))
    smiles['XLogP'] = smiles['mol'].apply(lambda x: Descriptors.MolLogP(x))
    smiles['NumValenceElectrons'] = smiles['mol'].apply(lambda x: Descriptors.NumValenceElectrons(x))              
    smiles['NumAromaticRings'] = smiles['mol'].apply(lambda x: Descriptors.NumAromaticRings(x))        
    smiles['NumSaturatedRings'] = smiles['mol'].apply(lambda x: Descriptors.NumSaturatedRings(x))        
    smiles['NumAliphaticRings'] = smiles['mol'].apply(lambda x: Descriptors.NumAliphaticRings(x))
    smiles['RingCount'] = smiles['mol'].apply(lambda x: Descriptors.RingCount(x))        
    smiles['TPSA'] = smiles['mol'].apply(lambda x: Descriptors.TPSA(x))
    smiles['mol'] = smiles['mol'].apply(lambda x: Chem.AddHs(x))
    smiles['num_of_atoms'] = smiles['mol'].apply(lambda x: x.GetNumAtoms())
    smiles['num_of_heavy_atoms'] = smiles['mol'].apply(lambda x: x.GetNumHeavyAtoms())
    c_patt = Chem.MolFromSmiles('C')
    def number_of_atoms(atom_list, data):
        for i in atom_list:
            data['num_of_{}_atoms'.format(i)] = data['mol'].apply(lambda x: len(x.GetSubstructMatches(Chem.MolFromSmiles(i))))

    number_of_atoms(['C','O', 'N', 'Cl'], smiles)
    return 0


"""
data = []
take_name = input("Drug name: \n")
take_name = take_name.replace(" ", "%20")
base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"+take_name+"/property/HBondDonorCount,HBondAcceptorCount,Charge,CanonicalSMILES,RotatableBondCount/JSON"
r = requests.get(base_url)
code = r.status_code
if code == 200:
    j = json.loads(r.content)
    data.append(j['PropertyTable']['Properties'])
else:
    print("Drug name not available in pubchem")
predict_data = pd.json_normalize(data[0])
mol_data(predict_data)"""



"""Make sure the main_data.csv file is in the same directory as the python file"""
mypath = os.path.dirname(os.path.abspath(__file__))
full_path = mypath+'\\main_data.csv'
model_data = pd.read_csv(full_path)
mol_data(model_data)



from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
X=model_data[['XLogP', 'Charge', 'HBondDonorCount', 'HBondAcceptorCount', 'RotatableBondCount', 'mol_w', 'num_of_atoms', 'num_of_C_atoms', 'num_of_O_atoms', 'num_of_N_atoms', 'num_of_Cl_atoms', 'NumValenceElectrons', 'TPSA']]
y=model_data['Phase']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3)

clf=RandomForestClassifier(n_estimators=100)
clf.fit(X_train,y_train)

Pkl_Filename = "RF_model.pkl"  

with open(Pkl_Filename, 'wb') as file:  
    pickle.dump(clf, file)

"""
y_pred=clf.predict(X_test)

print("Accuracy:",metrics.accuracy_score(y_test, y_pred))



prediction = clf.predict(predict_data[['XLogP', 'Charge', 'HBondDonorCount', 'HBondAcceptorCount', 'RotatableBondCount', 'mol_w', 'num_of_atoms', 'num_of_C_atoms', 'num_of_O_atoms', 'num_of_N_atoms', 'num_of_Cl_atoms', 'NumValenceElectrons', 'TPSA']])
predict_data['Phase'] = [int(prediction)]
predict_data['Drug'] = [take_name]
model_data = model_data.append(predict_data)
model_data.to_csv(full_path, index=False)


print("The given molecule is likely to enter Phase: ",int(prediction))




k=input("press enter to exit")"""


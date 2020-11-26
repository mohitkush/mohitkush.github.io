from flask import Flask, request, jsonify, render_template, url_for
import pandas as pd
import requests
import os
import pickle
from rdkit_model import mol_data
import json
import numpy as np


app = Flask(__name__)
model = pickle.load(open("RF_model.pkl", 'rb'))

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/predict', methods=['POST'])
def predict():
    data1 = []
    take_name = request.form
    drug_name = take_name['take_name']
    drug_name = str(drug_name.replace(" ", "%20"))
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"+str(drug_name)+"/property/HBondDonorCount,HBondAcceptorCount,Charge,CanonicalSMILES,RotatableBondCount/JSON"
    r = requests.get(base_url)
    code = r.status_code
    if code == 200:
        j = json.loads(r.content)
        data1.append(j['PropertyTable']['Properties'])
        print(data1)
    else:
        print("Drug name not available in pubchem")
    predict_data = pd.json_normalize(data1[0])
    mol_data(predict_data)
    prediction = model.predict(predict_data[['XLogP', 'Charge', 'HBondDonorCount', 'HBondAcceptorCount', 'RotatableBondCount', 'mol_w', 'num_of_atoms', 'num_of_C_atoms', 'num_of_O_atoms', 'num_of_N_atoms', 'num_of_Cl_atoms', 'NumValenceElectrons', 'TPSA']])
    return render_template('index.html' ,prediction_text="The given molecule is likely to enter Phase {}" .format(int(prediction)))

if __name__ == "__main__":
    app.run(debug=True)

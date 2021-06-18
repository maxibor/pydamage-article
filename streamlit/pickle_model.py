import pathlib
from pypmml import Model
import os 
import pickle

dirname = os.path.dirname(__file__)

model_path = dirname + "/../models/pydamage_glm_model.pmml"
print(model_path)
model = Model.load(model_path)
print(model.inputNames)
print(model.predict({'actual_cov':44.0, 'damage':0.01, 'contiglength':800}))
with open(dirname + "/../models/pydamage_glm_model_pmml.pickle", 'wb') as m:
    pickle.dump(model, m)
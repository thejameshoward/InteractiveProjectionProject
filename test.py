from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

import base64
from io import BytesIO
import PIL
from PIL import PngImagePlugin, Image

smiles = 'CC1=CC(=CC=C1)C(=O)CC(=O)C(F)(F)F'

mol = Chem.MolFromSmiles(smiles)

drawer = rdMolDraw2D.MolDraw2DCairo(500,180,200,180)
drawer.DrawMolecule(mol)
drawer.FinishDrawing()
im = BytesIO(drawer.GetDrawingText())
im = PngImagePlugin.PngImageFile(im)
im = Draw.MolsToGridImage([mol], subImgSize=(250,250))
buff = BytesIO()
im.save(buff, format='png')
encoded = base64.b64encode(buff.getvalue()).decode("utf-8")

old = 'data:image/jpeg;base64, ' + encoded

buff = BytesIO()
bio.save(buff, format='png')
encoded = base64.b64encode(buff.getvalue()).decode("utf-8")


new = 'data:image/jpeg;base64, ' + encoded



print(new)